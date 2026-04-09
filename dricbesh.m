function [xi] = dricbesh(nu,k,z)
%DRICBESJ 第三类球贝塞尔函数的导数（运用了递推公式来计算）
%   Detailed explanation goes here
if nargin == 2
    z = k;
    k = 1;
end
xi = 0.5*sqrt(pi/2/z)*besselh(nu+0.5,k,z)...
    + sqrt(pi*z/2)*dbesselh(nu+0.5,k,z);
end

