function [theta phi,r] = xcart2sph(x,y,z)
%XCART2SPH 直角坐标转换为球坐标
%   Detailed explanation goes here
theta = pi/2 - atan(z./sqrt(x.^2+y.^2));
phi = atan2(y,x);
r = sqrt(x.^2+y.^2+z.^2);
end

