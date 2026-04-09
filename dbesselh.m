function [h] = dbesselh(nu,k,z)
%DBESSELH 第三类贝塞尔函数的导数（运用了递推公式来计算）
%   Detailed explanation goes here
if nargin == 2
    z = k;
    k = 1;
end
h = 0.5*(besselh(nu-1,k,z) - besselh(nu+1,k,z));
end

