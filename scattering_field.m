clear all
close all
clc
parallel.gpu.enableCUDAForwardCompatibility(true)
%参数设置
%迭代求解器中的收敛阈值
CT=10^(-10);
%GPU选项
GPU=1;
%光波长
Lambda=532e-9;
%颗粒折射率：实部和虚部
np = 1.5;
%环境折射率
nb=1;
%半径
r_eff=1000e-9;
%分割
gird = 25;
%粒子形状
Np_shape= 'spherical';
% Np_shape='ellipsoid'; 
% Np_shape= 'rod';
% Np_shape= 'rec_block';
%粒子结构：单粒子还是双粒子
Structure="monomeric";
% Structure="dimeric";
%电场偏振方向和传播方向
E01=[0 0 1];%偏振方向为z
H01=[0 1 0]/(120*pi);
K01=[1 0 0];%传播方向为x正方向
%光场的形状
% IB="gaussian"; 
IB="plane wave" ;
% IB="laguerre-gaussian";
plane_2D="yz_x=2N";%前向散射
% plane_2D="yz_x=-2N";%后向散射
%侧向散射
% plane_2D="xz_y=2N";
% plane_2D="xz_y=2N";
% plane_2D="xy_z=2N";
% plane_2D="xy_z=-2N";
%坐标平面的场，观察粒子的电场增强现象
% plane_2D="xy" ;
% plane_2D="xz";
% plane_2D="yz";
%计算散射场
[Ex_out,Ey_out,Ez_out,Hx_out,Hy_out,Hz_out] = scattering(r_eff...
    ,gird,Lambda,Np_shape,Structure,IB,E01,H01,K01,CT,GPU,np,nb,plane_2D);
Sx = real((Ey_out.*conj(Hz_out)-Ez_out.*conj(Hy_out)));
Sy = real((Ez_out.*conj(Hx_out)-Ex_out.*conj(Hz_out)));
Sz = real((Ex_out.*conj(Hy_out)-Ey_out.*conj(Hx_out)));
%画图
figure;
subplot(331)
imagesc(abs(Ex_out))
title('Ex')
axis off
axis equal; % 保持纵横比
axis tight;
subplot(332)
imagesc(abs(Ey_out))
title('Ey')
axis off
axis equal; % 保持纵横比
axis tight;
subplot(333)
imagesc(abs(Ez_out))
title('Ez')
axis off
axis equal; % 保持纵横比
axis tight;
subplot(334)
imagesc(abs(Hx_out))
title('Hx')
axis off
axis equal; % 保持纵横比
axis tight;
subplot(335)
imagesc(abs(Hy_out))
title('Hy')
axis off
axis equal; % 保持纵横比
axis tight;
subplot(336)
imagesc(abs(Hz_out))
title('Hz')
axis off
axis equal; % 保持纵横比
axis tight;
subplot(337)
imagesc(Sx)
title('Sx')
axis off
axis equal; % 保持纵横比
axis tight;
subplot(338)
imagesc(Sy)
title('Sy')
axis off
axis equal; % 保持纵横比
axis tight;
subplot(339)
imagesc(Sz)
title('Sz')
axis off
axis equal; % 保持纵横比
axis tight;
% E_DDA = sqrt(Ex_out.*conj(Ex_out)+Ey_out.*conj(Ey_out)+Ez_out.*conj(Ez_out));
% H_DDA = sqrt(Hx_out.*conj(Hx_out)+Hy_out.*conj(Hy_out)+Hz_out.*conj(Hz_out));
% S_DDA = sqrt(Sx.*conj(Sx)+Sy.*conj(Sy)+Sz.*conj(Sz));