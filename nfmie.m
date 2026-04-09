function [E ,H] = nfmie(an ,bn,xc,yc,zc,r,nm,lambda)
%NFMIE 求给定展开系数的近场解E，H
%   Detailed explanation goes here
%ns为粒子折射率
%nm为背景折射率
%输入an，bn，xyz坐标，粒子半径，背景折射率，波长
k = 2*pi/lambda*nm;%背景介质中的波矢
M = numel(an);%an矩阵元素个数
n = (1:M)';%生成一个序列方便矩阵运算

E0 = 1;%入射波设置为均匀平面波
En = 1j.^n.*E0.*(2*n+1)./(n.*(n+1));%生成了一个序列

c0 = 299792458;
mue0 = 4*pi*1e-7;
omega = 2*pi/lambda*c0;%频率

E = zeros([numel(xc),3]);
H = zeros([numel(xc),3]);
R = zeros(3,3);
%有点C语言痕迹，
for ic = 1:numel(xc)
    [theta,phi,rad] = xcart2sph(xc(ic),yc(ic),zc(ic));%坐标变换
    rho = k*rad;
    [pin,taun] = angdepfun_mie(theta,n);
    h = sbesselh(n,rho);
    dxi = dricbesh(n,rho);
    %R为旋转矩阵
    R(1,1) = sin(theta)*cos(phi);
    R(1,2) = cos(theta)*cos(phi);
    R(1,3) = -sin(phi);
    R(2,1) = sin(theta)*sin(phi);
    R(2,2) = cos(theta)*sin(phi);
    R(2,3) = cos(phi);
    R(3,1) = cos(theta);
    R(3,2) = -sin(theta);
    R(3,3) = 0;
    
    M_oln = zeros(M,3);
    M_eln = zeros(M,3);
    N_oln = zeros(M,3);
    N_eln = zeros(M,3);
    
    Edum = zeros(M,3);
    Hdum = zeros(M,3);
    
    if rad > r
        M_oln(:,2) = cos (phi)*pin.* h;
        M_oln(:,3) = -sin (phi)*taun.*h ;
        M_eln(:,2) = -sin (phi)*pin.*h ;
        M_eln(:,3) = -cos (phi)*taun.*h ;
        N_oln( : , 1 ) = sin ( phi )*n.*(n+1)*sin(theta).*pin.* h/rho ;
        N_oln( : , 2 ) = sin ( phi ).*taun.*dxi/rho ;
        N_oln( : , 3 ) = cos ( phi ).*pin.*dxi/rho ;
        N_eln( : , 1 ) = cos ( phi )*n.*(n+1)*sin(theta).*pin.*h/rho ;
        N_eln( : , 2 ) = cos ( phi ).*taun.*dxi/rho ;
        N_eln( : , 3 ) = -sin ( phi ).*pin.*dxi/rho ;
        for d = 1:3
            Edum ( : , d ) = En .*( 1.j*an.'.* N_eln ( : , d ) - bn .'.* M_oln ( : , d ) ) ;
            Hdum ( : , d ) = En .*( 1.j*bn.'.* N_oln ( : , d ) + an .'.* M_eln ( : , d ) ) ; 
        end
    end
        E(ic,:) = R*sum(Edum).';
        H(ic,:) = k/omega/mue0*R*sum(Hdum).';
end
E = squeeze(reshape(E,[size(xc),3]));
H = squeeze(reshape(H,[size(xc),3]));
end

