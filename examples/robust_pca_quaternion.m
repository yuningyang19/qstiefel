function [U, cost] = robust_pca_quaternion(X, d)
% The following code is adapted from robust_pca.m in Manopt (available at
% https://www.manopt.org).    
% Author: Ying Wang, April, 2025.
% 
% Computes a robust version of PCA (principal component analysis) on
% quaternion data. 
% 
% function [U, cost] = robust_pca_quaternion(X, d)
%
% Given a quaternion matrix X of size p by n, such that each column
% represents a point in Q^p, this computes U: an orthonormal basis of size
% p by d such that the column space of U captures the points X as well as
% possible. More precisely, the function attempts to compute U as the
% minimizer over the quaternion Stiefel manifold of:
%
%  f(U) = (1/n) Sum_{i = 1:n} dist(X(:, i), the space spanned by U)
%       = (1/n) Sum_{i = 1:n} || U*U'*X(:, i) - X(:, i) ||
%
% The output cost represents the average distance achieved with the
% returned U. Notice that norms are not squared, for robustness. 
%
% In practice, because this function is nonsmooth, it is smoothed with a
% pseudo-Huber loss function of parameter epsilon (noted e for short), and
% the smoothing parameter is iteratively reduced (with warm starts):
%
%   f_e(U) = (1/n) Sum_{i = 1:n} l_e(|| U*U'*X(:, i) - X(:, i) ||)
%
%   with l_e(x) = sqrt(x^2 + e^2) - e (for e = 0, this is absolute value).
%
% The intermediate optimization of the smooth cost over the quaternion
% Stiefel manifold is performed using the Manopt toolbox and the file
% stiefelquaternionfactory.m in this package. 
%
% Ideally, the non-outlier data should be centered. If not, this
% pre-processing centers all the data, but bear in mind that outliers will
% shift the center of mass too.
% X = X - repmat(mean(X, 2), [1, size(X, 2)]);
%
% There are no guarantees that this code will return the optimal U.

    % Prepare a Manopt problem structure for optimization of the given
    % cost (defined below) over the quaternion Stiefel manifold.
    [p, n] = size(X);
    manifold = stiefelquaternionfactory(p, d);
    problem.M = manifold;
    problem.cost = @robustpca_cost;
    problem.egrad = @robustpca_gradient;	
    
	% Iteratively reduce the smoothing constant epsilon and optimize
	% the cost function over quaterion Stiefel manifold.
    epsilon = 1;
	n_iterations = 6;
	reduction = .5;
	options.verbosity = 0; % Change this number for more or less output
    warning('off', 'manopt:getHessian:approx');    
    U = quaternion(qr_unique(randn(p,d)));
    
    for iter = 1 : n_iterations
        fprintf('iteration %d\n',iter);
        U = trustregions(problem, U, options);
        epsilon = epsilon * reduction;
    end
    warning('on', 'manopt:getHessian:approx');    
    
	% Return the cost as the actual sum of distances, not smoothed.
	epsilon = 0;
	cost = robustpca_cost(U);          
    
    % Smoothed cost
    function value = robustpca_cost(U)

        vecs = U*(U'*X) - X;
        sqnrms = sum(abs(vecs).^2, 1);
        vals = sqrt(sqnrms + epsilon^2) - epsilon;
        value = mean(vals);

    end

    % Euclidean gradient of the smoothed cost (it will be transformed into
    % the Riemannian gradient automatically by Manopt).
    function G = robustpca_gradient(U)

        UtX = U'*X;
        vecs = U*UtX-X;
        sqnrms = sum(abs(vecs).^2, 1);
        G = zeros(p, d);
        for i = 1:n
            G = G + (1/sqrt(sqnrms(i) + epsilon^2)) * vecs(:,i) * UtX(:,i)';
        end
        G = G/n;
    end
end