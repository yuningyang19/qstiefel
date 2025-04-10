function [U, S, V, info] = truncated_svd_quaternion(A, p)
% The following code is adapted from truncated_svd.m in Manopt (available
% at https://www.manopt.org).   
% Author: Ying Wang, April, 2025.
% 
% Returns a rank p truncated SVD of a quaternion matrix A. 
%
% function [U, S, V, info] = truncated_svd(A, p)
%
% Input: A quaternion matrix A of size mxn and an integer p <= min(m, n).
% Output: An orthonormal quaternion matrix U of size mxp, an orthonormal
%         quaternion matrix Y of size nxp and a diagonal real matrix S of
%         size pxp with nonnegative and decreasing diagonal entries such
%         that USV.' is the best rank p approximation of A according to the
%         Frobenius norm. This function produces an output akin to svds.
% 
% The optimization is performed on the quaternion Stiefel manifold and uses
% the code stiefelquaternionfactory.m provided in this package. 
% 
% The decomposition is obtained by maximizing
%   f(U, V) = .5*norm(U'*A*V, 'fro')^2
% where U, V are orthonormal. It is easy to show that maximizing f is
% equivalent to minimizing g with 
%   g(U, V) = min_S norm(U*S*V' - A, 'fro')^2,
% which confirms that we are going for a best low-rank approximation of A.
% 
% The code can be modified to accept a function handle for A(x) = A*x
% instead of a matrix A, which is often useful. This would further require
% a function handle At for the Hermitian transpose of A, such that At(x) = A'*x.
    
    % Generate some random data to test the function if none is given.
    if ~exist('A', 'var') || isempty(A)
        A = randq(42, 60); % randq is in QTFM.
    end
    if ~exist('p', 'var') || isempty(p)
        p = 5;
    end
    
    % Retrieve the size of the problem and make sure the requested
    % approximation rank is at most the maximum possible rank.
    [m, n] = size(A);
    assert(p <= min(m, n), 'p must be smaller than the smallest dimension of A.');
    
    % Define the cost and its derivatives on the quaternion Stiefel manifold
    tuple.U = stiefelquaternionfactory(m, p);
    tuple.V = stiefelquaternionfactory(n, p);
    M = productmanifold(tuple);
    
    % Define a problem structure, specifying the manifold M, the cost
    % function and its derivatives.
    problem.M = M;
    problem.cost  = @cost;
    problem.egrad = @egrad;
    problem.ehess = @ehess;
    
    % The functions below make many redundant computations. This
    % performance hit can be alleviated by using the caching system.
    
    % Cost function
    function f = cost(X)
        U = X.U;
        V = X.V;
        f = -.5*normQf(U'*A*V)^2; % normQf is in Qlab.
    end
    % Euclidean gradient of the cost function
    function g = egrad(X)
        U = X.U;
        V = X.V;
        AV = A*V;
        AtU = A'*U;
        g.U = -AV*(AV'*U);
        g.V = -AtU*(AtU'*V);
    end
    % Euclidean Hessian of the cost function
    function h = ehess(X, H)
        U = X.U;
        V = X.V;
        Udot = H.U;
        Vdot = H.V;
        AV = A*V;
        AtU = A'*U;
        AVdot = A*Vdot;
        AtUdot = A'*Udot;
        h.U = -(AVdot*AV'*U + AV*AVdot'*U + AV*AV'*Udot);
        h.V = -(AtUdot*AtU'*V + AtU*AtUdot'*V + AtU*AtU'*Vdot);
    end
    
    % Execute some checks on the derivatives for early debugging.
    % These things can be commented out of course.
    % checkgradient(problem);
    % pause;
    % checkhessian(problem);
    % pause;
    
    % Issue a call to a solver. A random initial guess will be chosen and
    % default options are selected. 
    options.Delta_bar = 4*sqrt(2*p);
    options.tolgradnorm = 1e-9;
    options.verbosity = 0;
    [X, Xcost, info] = trustregions(problem, [], options); %#ok<ASGLU>
    U = X.U;
    V = X.V;
    
    % Finish the job by rotating U and V such that the middle matrix S can
    % be diagonal with nonnegative, decreasing entries. This requires a
    % small svd of size pxp.
    Spp = U'*A*V;
    [Upp, Spp, Vpp] = svd(Spp);
    U = U*Upp;
    S = Spp;
    V = V*Vpp;    
end