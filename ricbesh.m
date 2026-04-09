function [xi] = ricbesh(nu,k,z)
%RICHBESH 第三类贝塞尔函数
%   Detailed explanation goes here
if nargin == 2
    z = k;
    k = 1;
end
xi = sqrt(pi*z/2)*besselh(nu+0.5,k,z);
end

