function [X, A] = PCA_stochastic_quaternion(A, k)
% The following code is adapted from PCA_stochastic.m in Manopt (available
% at https://www.manopt.org).   
% Author: Ying Wang, April, 2025.
% 
% Example of stochastic gradient algorithm for a quaternion PCA problem.
% 
% QPCA (quaternion principal component analysis) on a quaternion matrix A
% of size nxd consists in solving
% 
%   minimize_X  f(X) = -.5*norm(A*X, 'fro')^2 / n,
% 
% where X is a quaternion matrix of dimension dxk with orthonormal columns.
% This is equivalent to finding k dominant singular vectors of A, or k top
% eigenvectors of A'*A. If n is large, this computation can be expensive.
% Thus, stochastic gradient algorithms take the point of view that f(X) is
% a sum of many (n) terms: each term involves only one of the n rows of A.
%
% To make progress, it may be sufficient to optimize with respect to a
% subset of the terms at each iteration. This way, each individual
% iteration can be very cheap. In particular, individual operations have
% cost independent of n, because f or its gradient need never be evaluated
% completely (or at all in the case of f.)
%
% Stochastic gradient algorithms (this implementation in particular) are
% sensitive to proper parameter tuning. See in code.

    % If none is given, generate a random data set: n samples in Q^d
    if ~exist('A', 'var') || isempty(A)
        d = 1000;
        n = 100000;
        fprintf('Generating data...');
        A = randq(n, d)*quaternion(diag([[15 10 5], ones(1, d-3)])); % randq is in QTFM
        fprintf(' done (size: %d x %d).\n', size(A));
    else
        [n, d] = size(A);
    end

    % Pick a number of component to compute
    if ~exist('k', 'var') || isempty(k)
        k = 3;
    end
    
    % We are looking for k orthonormal vectors in Q^d: quaternion Stiefel manifold.
    problem.M = stiefelquaternionfactory(d, k);
    
    % The cost function to minimize is a sum of n terms. This parameter
    % must be set for stochastic algorithms.
    problem.ncostterms = n;
    
    % We do not need to specify how to compute the value of the cost
    % function (stochastic algorithms never use this). All we need is to
    % specify how to compute the gradient of the cost function, where the
    % sum is restricted to a subset of the terms (a sample). Notice that we
    % specify a partial Euclidean gradient (hence the 'e' in partialegrad).
    % This way, Manopt will automatically convert the Euclidean vector into
    % a proper Riemannian partial gradient, in the tangent space at X.
    % In particular, if sample = 1:n, then the partial gradient corresponds
    % to the actual (complete) gradient.
    problem.partialegrad = @partialegrad;
    function G = partialegrad(X, sample)
        
        % X is a quaternion orthonormal matrix of size dxk
        % sample is a vector of indices between 1 and n
        % Extract a subset of the dataset
        Asample = A(sample, :);
        
        % Compute the gradient of f restricted to that sample
        G = -Asample'*(Asample*X);
        G = G / n;
        
    end

    % If one wants to use checkgradient to verify one's work, then it is
    % necessary to specify the cost function as well, as below.
    % problem.cost = @(X) -.5*normQf(A*X)^2 / n; % normQf is in Qlab
    % checkgradient(problem); pause;

    % To have the solver record statistics every x iterations, set
    % options.checkperiod to x. This will record simple quantities which
    % are almost free to compute (namely, elapsed time and step size of the
    % last step.) To record more sophisticated quantities, you can use
    % options.statsfun as usual. Time spent computing these statistics is
    % not counted in times reported in the info structure returned by the
    % solver.
    options.checkperiod = 10;
    options.statsfun = statsfunhelper('metric', @(X) normQf(A*X));
    
    % Set the parameters for the solver: stochastic gradient algorithms
    % tend to be quite sensitive to proper tuning, especially regarding
    % step size selection. See the solver's documentation (in Manopt) for
    % details. 
    options.maxiter = 1000;
    options.batchsize = 10;
    % options.stepsize_type = 'decay';
    options.stepsize_init = 1e2;
    options.stepsize_lambda = 1e-3;
    options.verbosity = 1;
    
    % Run the solver
    [X, info] = stochasticgradient(problem, [], options);    
    
    % Plot the special metric recorded by options.statsfun
    plot([info.iter], [info.metric], '.-');
    xlabel('Iteration #');
    ylabel('Frobenius norm of A*X');
    title('Convergence of stochasticgradient on stiefelquaternionfactory for quaternion PCA');
    
    % Add to that plot a reference: the globally optimal value attained if
    % the true dominant eigenvectors are computed by calling the code
    % dominant_invariant_subspace_quaternion.m in the qstiefel package. 
    fprintf('Running dominant_invariant_subspace_quaternion for comparison... ');
    t = tic();
    B = A'*A; B = .5*(B+B');
    [V, ~, ~] = dominant_invariant_subspace_quaternion(B, k);
    fprintf('done: %g [s]\n', toc(t));
    bound = normQf(A*V); arraymetric = [info.metric];
    fprintf('Difference in objective value: %d\n', bound-arraymetric(end));
    hold all;    
    plot([info.iter], bound*ones(size([info.iter])), '--');
    hold off;    
    legend('Algorithm', 'EVD bound', 'Location', 'SouthEast');        
end