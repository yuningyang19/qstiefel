function [X, iter] = svtQAlg(M, Omega, tau, delta, max_iter, tol)
% singular value thresholding algorithm for quaternion matrix completion,
% a simple implementation for testing

% M: observed matrix with 0 in unobserved positions
% Omega: binary matrix indicating observed positions (1 for observed, 0 for unobserved)
% tau: thresholding parameter
% delta: step size
% max_iter: maximum number of iterations
% tol: tolerance for convergence

% Initialize variables
[m, n] = size(M);
Y = M;
 
    M_Omega = M .* Omega;
% Iterate
for iter = 1:max_iter
    % Singular value thresholding
    [X] = svtQ(Y, tau);
 

    % Update Y
    X_Omega = X .* Omega;

    Y = Y + delta * (M_Omega - X_Omega);

    % Check convergence
    error = normQf(M_Omega - X_Omega) / normQf(M_Omega);
    fprintf('iter: %d, error = %.4f\n', iter, error)
    if error < tol
        break;
    end
end
end



