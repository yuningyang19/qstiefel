function R = sprandnQ(m,n,dens)
%SPRANDNQ 此处显示有关此函数的摘要
%   此处显示详细说明

w = full(sprandn(m,n,dens)); 
x = full(sprandn(m,n,dens)); 
y = full(sprandn(m,n,dens)); 
z = full(sprandn(m,n,dens)); 
R = quaternion(w,x,y,z);
end

