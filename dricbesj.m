function [S] = dricbesj(nu,z)
%DRICBESJ 第一类球贝塞尔函数的导数（运用了递推公式来计算）
%   Detailed explanation goes here
S = 0.5*sqrt(pi/2/z)*besselj(nu+0.5,z) + sqrt(pi*z/2)*dbesselj(nu+0.5,z);
end

