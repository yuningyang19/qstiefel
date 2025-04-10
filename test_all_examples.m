% Test of optimization over quaternion Stiefel manifold
clear;
eigs = true; svds = true; sparsepca = true; robustpca = true; stochasticpca = true;
%% Dominant Eigenvalue Decomposition
if eigs
    n = 1000; p = 5;
    A = randq(n); % one can also use: A = quaternion(rand(n),rand(n),rand(n),rand(n));
    A = 0.5*(A+A');
    tic; [X,~,~] = dominant_invariant_subspace_quaternion(A, p); time = toc; rienorm = rienorm_eigs(A,X);
    fprintf('\ndominant eigenvalue decomposition\n\nqstiefel (calling Manopt):\ntime: %.2fs\nRiemannian norm: =%d\n',time,rienorm);
    tic; [X_qtfm,S_qtfm] = eig(A); [~,ind] = sort(diag(S_qtfm),'descend'); X_qtfm = X_qtfm(:,ind); X_qtfm = X_qtfm(:,1:p); time_qtfm = toc; rienorm_qtfm = rienorm_eigs(A,X_qtfm);
    fprintf('QTFM:\ntime: %.2fs\nRiemannian norm: %d\n',time_qtfm,rienorm_qtfm);
end
%% Dominant Singular Value Decomposition
if svds
    m = 1000; n = 1000; p = 5;
    A = randq(m,n);
    tic; [U, S, V, ~] = truncated_svd_quaternion(A, p); time = toc; rienorm = rienorm_svds(A,U,V);
    fprintf('\ndominant singular value decomposition\n\nqstiefel (calling Manopt):\ntime: %.2fs\nRiemannian norm: =%d\n',time,rienorm);
    tic; [U_qtfm,~,V_qtfm] = svd(A); U_qtfm = U_qtfm(:,1:p); V_qtfm = V_qtfm(:,1:p); time_qtfm = toc; rienorm_qtfm = rienorm_svds(A,U_qtfm,V_qtfm);
    fprintf('QTFM:\ntime: %.2fs\nRiemannian norm: %d\n',time_qtfm,rienorm_qtfm);
end
%% Sparse Quaternion PCA
if sparsepca
    m = 100; n = 50; p = 5;
    A = randq(m,n);
    gamma = 1;
    fprintf('\nSparse Quaternion PCA begins:\n');
    tic; [Z, P, ~] = sparse_pca_quaternion(A, p, gamma); time = toc;
    fprintf('\nSparse Quaternion PCA completed in %.2fs.\n',time);
end
%% Robust Quaternion PCA
if robustpca
    m = 50; n = 100; p = 5;
    A = randq(m,n);
    fprintf('\nRobust Quaternion PCA begins:\n\n');
    tic; [U, cost] = robust_pca_quaternion(A, p); time = toc;
    fprintf('\nRobust Quaternion PCA completed in %.2fs.\n',time);
end
%% Stochastic Quaternion PCA
if stochasticpca
    fprintf('\nStochastic Quaternion PCA begins:\n');
    [X, ~] = PCA_stochastic_quaternion();
end
%% Auxiliary Functions
function norm = rienorm_eigs(A,U)
G = -A*U;
gradU = G - U*herm(G'*U);
norm = normQf(gradU);
end

function norm = rienorm_svds(A,U,V)
AV = A*V;
UstarAV = U'*AV;
gradU = -AV + U*herm(UstarAV);
gradV = -A'*U+V*herm(UstarAV');
norm = normQf(gradU)^2 + normQf(gradV)^2; norm = sqrt(norm);
end

function Aherm = herm(A)
Aherm = 0.5*(A+A');
end