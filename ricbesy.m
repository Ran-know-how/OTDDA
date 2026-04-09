function [C] = ricbesy(nu,z)
%RICBESY 第二类贝塞尔函数
%   Detailed explanation goes here
S = -sqrt(pi*z/2)*bessely(nu+0.5,z);
end

