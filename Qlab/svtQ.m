function [X] = svtQ(A,tau,varargin)
%SVTQ the singular value thresholding operator for quaternion matrices
% varargin{1}: the qsvd method, default = @ModifiedcsvdQ77

    if ~isempty(varargin)
        qsvd = varargin{1};
    else
        qsvd = @ModifiedcsvdQ77;
    end
    
    %qsvd = @csvdQ;
    tic,
    [u,s,v] = qsvd(A);
    % [u,s,v] = ranQLoRMA(A,floor(size(A,1)/3));
    t = toc;
    stau = max(s-tau,0);
    X = u*stau*v';

    fprintf("Ours modified complex based svd:\n")
    fprintf("residual of ||b-u*s*v||: %.15e, time = %.2f\n", normQf(A-u*s*v'), t)
%fprintf("residual of singular values:  %.15e\n", norm(s3-s))
fprintf("cond of U and V: %.15e, %.15e\n\n", condQ(u), condQ(v))



end

