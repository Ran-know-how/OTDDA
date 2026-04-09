function [E,H] = calcmie_nf(r,ns,nm,lambda,xc,yc,zc)
%CALCMIE_NF 计算单球电磁散射近场
%ns为粒子折射率
%nm为背景折射率

k = 2*pi/lambda*nm;
x = k*r;
m = ns/nm;

[an,bn] = expcoeff_mie(x,m);

[E,H] = nfmie(an,bn,xc,yc,zc,r,nm,lambda);
end

