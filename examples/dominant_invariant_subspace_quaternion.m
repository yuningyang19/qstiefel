function [X, D, info] = dominant_invariant_subspace_quaternion(A, p)
% The following code is adapted from dominant_invariant_subspace.m in
% Manopt (available at https://www.manopt.org).  
% Author: Ying Wang, April, 2025.
% 
% Returns an orthonormal basis of the dominant invariant p-subspace of a
% quaternion Hermitian matrix A. 
%
% function [X, D, info] = dominant_invariant_subspace(A, p)
%
% Input: A quaternion, Hermitian matrix A of size nxn and an integer p <= n.
% Output: A quaternion, orthonormal matrix X of size nxp such that
%         re(trace(X'*A*X)) is maximized. That is, the columns of X form an
%         orthonormal basis of a dominant invariant subspace of dimension p
%         of A. We also perform a Rayleigh-Ritz procedure, so eigenvectors
%         associated with the largest eigenvalues of A are obtained. Sign
%         is important: 2 is deemed a larger eigenvalue than -5. D contains
%         the largest p eigenvalues.
%
% The optimization is performed on the quaternion Stiefel manifold and uses
% the code stiefelquaternionfactory.m provided in this package. 

    % Generate some random data to test the function
    if ~exist('A', 'var') || isempty(A)
        A = randq(128); % randq is in QTFM.
        A = (A+A')/2;
    end
    if ~exist('p', 'var') || isempty(p)
        p = 3;
    end
    
    % Make sure the input matrix is square and Hermitian
    n = size(A, 1);
    assert(size(A, 2) == n, 'A must be square.');
    assert(normQf(A-A') < n*eps, 'A must be Hermitian.'); % normQf is in Qlab.
	assert(p<=n, 'p must be smaller than n.');
    
    % Define the cost and its derivatives on the quaternion Stiefel manifold
    St = stiefelquaternionfactory(n, p);
    problem.M = St;
    problem.cost = @(X)    -.5*(sum(dot(X,A*X)).w);
    problem.grad = @(X)    -St.egrad2rgrad(X, A*X);
    problem.hess = @(X, H) -St.ehess2rhess(X, A*X, A*H, H);    
    
    % Execute some checks on the derivatives for early debugging.
    % These can be commented out.
    % checkgradient(problem);
    % pause;
    % checkhessian(problem);
    % pause;
    
    % Issue a call to a solver. A random initial guess will be chosen and
    % default options are selected except for the ones we specify here.
    options.Delta_bar = 8*sqrt(p);
    options.tolgradnorm = 1e-9;
    options.verbosity = 0;
    [X, costX, info] = trustregions(problem, [], options); %#ok<ASGLU>

    % Rayleigh-Ritz procedure.
    S = X'*A*X; S = .5*(S+S');
    [U, S] = eig(S);
    X = X*U;
    D = S;    
end