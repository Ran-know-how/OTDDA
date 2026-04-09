function [S] = ricbesj(nu,z)
%RICBESJ 第一类贝塞尔函数
%   Detailed explanation goes here
%球贝塞尔函数
S = sqrt(pi*z/2)*besselj(nu+0.5,z);
end

