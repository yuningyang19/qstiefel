function R = sprandQ(m,n,density)
%SPRANDQ Summary of this function goes here
%   Detailed explanation goes here

    R = quaternion(full(sprand(m,n,density)),full(sprand(m,n,density)),full(sprand(m,n,density)),full(sprand(m,n,density)));

end

