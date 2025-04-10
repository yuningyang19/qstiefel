function M = stiefelquaternionfactory(n, p, k)
% The following code is adapted from stiefelfactory.m in Manopt (available
% at https://www.manopt.org). 
% The code relies on Quaternion toolbox for Matlab (QTFM, available at
% https://qtfm.sourceforge.io).
% The code calls several functions in Qlab (by Chao Chang and Yuning Yang).
% Author: Ying Wang, April, 2025.
% 
% Returns a manifold struct. to optimize over quaternion orthonormal matrices.
%
% function M = stiefelquaternionfactory(n, p)
% function M = stiefelquaternionfactory(n, p, k)
% 
% The quaternion Stiefel manifold is the set of quaternion orthonormal nxp
% matrices. If k is larger than 1, this is the Cartesian product of the
% quaternion Stiefel manifold taken k times. The metric is such that the
% manifold is a Riemannian submanifold of Q^nxp (where Q is the skew-field 
% of quaternions) equipped with the real-trace inner product (R-product).
%
% Points are represented as matrices X of size n x p x k (or n x p if k=1,
% which is the default) such that each n x p quaternion matrix is orthonormal,
% i.e., X'*X = eye(p) if k = 1, or X(:, :, i)' * X(:, :, i) = eye(p) for
% i = 1 : k if k > 1. Tangent vectors are represented as matrices the same
% size as points.
%
% The default retraction is the polar retraction (instead of the QR
% retraction). This is motivated by the observation that available
% algorithms for quaternion QR decomposition are not as efficient as their
% real/complex counterparts while the quaternion polar retraction can be
% almost as efficient as its real/complex counterparts since it is computed
% by solving quaternion linear equations, which can be reduced to solving
% real/complex linear equations. 
% 
% To use the QR retraction, run
%    M.retr = M.retr_qr;
% after creating M with this factory. This can be reverted with
%    M.retr = M.retr_polar.
% 
% To use the Cayley retraction, run
%    M.retr = M.retr_cayley;
% after creating M with this factory. This can be reverted with
%    M.retr = M.retr_polar.
%
% By default, k = 1. 

    assert(n >= p, 'The dimension n must be larger than the dimension p.');
    
    if ~exist('k', 'var') || isempty(k)
        k = 1;
    end

    array_type = 'double';      
    
    if k == 1
        M.name = @() sprintf('Quaternion Stiefel manifold St_Q(%d, %d)', n, p);
    elseif k > 1
        M.name = @() sprintf('Product Quaternion Stiefel manifold St_Q(%d, %d)^%d', n, p, k);
    else
        error('k must be an integer no less than 1.');
    end
    
    M.dim = @() k*(4*n*p - 2*p*p + p);
    
    M.inner = @(x, d1, d2) sum(dot(d1,d2)).w; % equivalent to: sum(d1(:)'*d2(:)).w;
    
    M.norm = @(x, d) normQf(d); % normQf is in Qlab.
    
%     M.typicaldist = @() sqrt(p*k);
    
    M.proj = @projection;
    function Up = projection(X, U)
        % The following code is not yet implementable since multiprod is not supported in QTFM.
        %         XtU = multiprod(multitransp(X), U);
        %         symXtU = multisym(XtU);
        %         Up = U - multiprod(X, symXtU);
        % The code above, when implementable, shoule be faster than the following equivalent code.
        Up = zeros(size(U)); Up = quaternion(Up); % convert Up to quaternion format.
        for i = 1 : k
            Xi = X(:, :, i);
            Ui = U(:, :, i);
            Up(:, :, i) = Ui - Xi*herm(Xi'*Ui);
        end
    end
    
    M.tangent = M.proj;
    
    M.tangent2ambient_is_identity = true;

    M.tangent2ambient = @(X, U) U;
    
    % For Riemannian submanifolds, converting a Euclidean gradient into a
    % Riemannian gradient amounts to an orthogonal projection.
    M.egrad2rgrad = M.proj;
    
    M.ehess2rhess = @ehess2rhess;
    function rhess = ehess2rhess(X, egrad, ehess, H)
        %         XtG = multiprod(multitransp(X), egrad);
        %         symXtG = multisym(XtG);
        %         HsymXtG = multiprod(H, symXtG);
        %         rhess = projection(X, ehess - HsymXtG);
        rhess = zeros(size(H)); rhess = quaternion(rhess);
        for i = 1 : k
            Xi = X(:,:,i);
            egradi = egrad(:,:,i);
            ehessi = ehess(:,:,i);
            Hi = H(:,:,i);
            XtGi = Xi'*egradi;
            hermXtGi = herm(XtGi);
            HsymXtGi = Hi*hermXtGi;
            rhess(:,:,i) = projection(Xi, ehessi - HsymXtGi);
        end
    end
    
    M.retr_qr = @retraction_qr;
    function Y = retraction_qr(X, U, t)
        Y = zeros(size(X)); Y = quaternion(Y);
        for i = 1:k
            Xi = X(:,:,i); Ui = U(:,:,i);
            if nargin < 3
                Yi = qr(Xi + Ui); Yi = Yi(:,1:p); Y(:,:,i) = Yi; % One can use other quaternion qr codes.
            else
                Yi = qr(Xi + t*Ui); Yi = Yi(:,1:p); Y(:,:,i) = Yi;
            end
        end
    end
    
    M.retr_polar = @retraction_polar;
    function Y = retraction_polar(X, U, t)
        if nargin < 3
            Y = X + U;
        else
            Y = X + t*U;
            U = t*U;
        end
        for kk = 1 : k
%             % The following code is not implementable yet since economic-size svd is not supported in QTFM yet.
%             [u, s, v] = svd(Y(:, :, kk), 'econ'); Y(:, :, kk) = u*v';
            b = Y(:, :, kk);
            Ukk = U(:, :, kk); 
            A = eye(p) + Ukk'*Ukk;
            A = sqrtmQ(A);
            Y(:, :, kk) = QLEQL(b,A); % This function is in Qlab and solves the linear equation Y*A=b.
        end
    end    

    M.retr_cayley = @retraction_cayley;
    function Y = retraction_cayley(X, U, t)
        % Y = (I-0.5A)^{-1}(I+0.5A)X where A = PUX'-XU'P with P=I-0.5XX'.
        % We use the following low-rank version.
        if nargin == 3
            U = t*U;
        end
        Y = quaternion(zeros(size(X)));
        for kk = 1 : k
            Xkk = X(:, :, kk);
            Ukk = U(:, :, kk);
            PU = Ukk-0.5*Xkk*(Xkk'*Ukk); Utemp = [PU,Xkk]; Vtemp = [Xkk,-PU];
            A = eye(2*p)-0.5*Vtemp'*Utemp;
            b = Vtemp'*Xkk;
            S = QLEQ(A,b); % This function is in Qlab and solves the linear equation A*S=b.
            Y(:, :, kk) = Xkk + Utemp*S;
        end
    end    
    
    % By default, we use the polar retraction.
    M.retr = M.retr_polar;

    % The polar retraction is second order, and tagged as such in retr2.
    M.retr2 = M.retr_polar;

    M.exp = @exponential;
    function Y = exponential(X, U, t)
        if nargin == 2
            tU = U;
        else
            tU = t*U;
        end
        Y = zeros(size(X), array_type); Y = quaternion(Y);
        I = eye(p, array_type); I = quaternion(I);
        Z = zeros(p, array_type); Z = quaternion(Z);
        for kk = 1 : k
            Xkk = X(:, :, kk);
            Ukk = tU(:, :, kk);
            Y(:, :, kk) = [Xkk Ukk] * ...
                         expm([Xkk'*Ukk , -Ukk'*Ukk ; I , Xkk'*Ukk]) * ...
                         [ expm(-Xkk'*Ukk) ; Z ];
        end        
    end
    
    M.rand = @randqst;
    function X = randqst()
        X = quaternion(qr_unique(randn(n, p, k, array_type)));
%         X = randq(n, p, k);
%         for kk = 1:k
%             [Xkk,~] = qr(X(:, :, kk)); Xkk = Xkk(:, 1:p); X(:, :, kk) = Xkk;            
%         end
    end 
    
    M.randvec = @randomvec; % random tangent vector
    function U = randomvec(X)
        U = projection(X, quaternion(randn(n, p, k, array_type)));
        U = U / normQf(U);
    end
    
    M.lincomb = @matrixlincomb;
    
    M.zerovec = @(x) quaternion(zeros(n, p, k, array_type));
    
    M.transp = @(x1, x2, d) projection(x2, d);
    
    M.vec = @(x, u_mat) u_mat(:);
    M.mat = @(x, u_vec) reshape(u_vec, [n, p, k]);
    M.vecmatareisometries = @() true;      
end

function A = herm(A), A = 0.5*(A+A'); end