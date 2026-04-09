function [C] = dricbesy(nu,z)
%DRICBESJ 第二类球贝塞尔函数的导数（运用了递推公式来计算）
%   Detailed explanation goes here
C = -0.5*sqrt(pi/2/z)*bessely(nu+0.5,z) + sqrt(pi*z/2)*dbessely(nu+0.5,z)
end

