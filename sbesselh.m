function [h] = sbesselh(nu,k,z)
%SBESSELH 计算第三类球贝塞尔函数
%   Detailed explanation goes here
if nargin == 2
    z = k;
    k = 1;
end
h = ricbesh(nu,k,z)./z;
end

