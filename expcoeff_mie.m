function [an , bn] = expcoeff_mie(x,m)
%与expcoeff_mie相比多输出一阶an，bn，用在计算g中
M = ceil(x + 4*(x^(1/3))+2)+1;%截断系数,与尺寸参数有关
n = 1:M;
Sx = ricbesj(n,x);%第一类球贝塞尔函数
dSx = dricbesj(n,x);%第一类球贝塞尔函数的倒数
 
Smx = ricbesj(n,m*x);%第一类球贝塞尔函数，mx处
dSmx = dricbesj(n,m*x);%第一类球贝塞尔函数的导数，mx处
 
xix = ricbesh(n,x);%第三类球贝塞尔函数
dxix = dricbesh(n,x);%第三类球贝塞尔函数的导数
 
an = (m*Smx.*dSx - Sx.*dSmx)./(m*Smx.*dxix - xix.*dSmx);%an公式
bn = (Smx.*dSx - m*Sx.*dSmx)./(Smx.*dxix - m*xix.*dSmx);%bn公式
end

