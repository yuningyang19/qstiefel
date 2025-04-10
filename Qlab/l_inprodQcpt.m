function [p] = l_inprodQcpt(a,b)
%L_INPRODCPT left inner product a.'*conj(b) between two quaternion matrices (vectors)
% a and b represented compactly as a = [a0 a1], b = [b0 b1] with a = a0+a1j
% and b = b0 + b1j
if mod(size(a,2),2) == 1 || mod(size(b,2),2) == 1   
    error('not a compact representation of a quaternion matrix')
end
if size(a,1) ~= size(b,1) || size(a,2) ~= size(b,2)
    error('size not matched')
end

p = mtimesQcpt(transQcpt(mat2vecQcpt(a)),conjQcpt(mat2vecQcpt(b)));
end

