clear
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
np = 1.5+0.1i;
%环境折射率
nb=1;
%半径
r_eff=500e-9;
%分割
gird = 10;
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
IB="gaussian";
% IB="plane wave" ;
% IB="laguerre-gaussian";
if IB=="plane wave" 
    z0=0;              % Focus point of the Gaussian beam
    Waist_r=100;       % ratio of waist raduis of beam to wavelength
    p = 0;
    l= 0;
elseif IB=="gaussian"
    z0=0;              % Focus point of the Gaussian beam
    NA = 0.3;
    n_object = 1;
    Waist_r = n_object/(pi*NA);
    p = 0;
    l= 0;
elseif IB=="laguerre-gaussian"
    z0=0;              % Focus point of the Gaussian beam
    NA = 0.1;
    n_object = 1;
    Waist_r = n_object/(pi*NA);     
    p = 3;
    l= 0;
end 
%力的遍历方向
direction = "y";
%遍历的范围
l_scan = -3*r_eff : r_eff/2 : 3*r_eff;              % 粒子中心z坐标遍历数组
[F_line,T_line] = force_line(r_eff,gird,Lambda,Np_shape...
    ,Structure,IB,E01,H01,K01,CT,GPU,np,nb,direction,l_scan,z0,p,l,Waist_r);
% 6. 结果格式转换（GPU→CPU，便于绘图）
if GPU == 1
    F_line = gather(F_line);
    T_line = gather(T_line);
    l_scan = gather(l_scan);
end
figure('Position', [100, 100, 800, 400]);
plot(l_scan, F_line);
legend Fx Fy Fz
xlabel('粒子位置（m）');
ylabel('力（N）');
title('光镊纵向力（国际单位制）');
grid on;
figure('Position', [100, 100, 800, 400]);
plot(l_scan, T_line);
legend Tx Ty Tz
xlabel('粒子位置（m）');
ylabel('力矩（N*m）');
title('光镊的力矩（国际单位制）');
grid on;
% set(gca, 'YScale', 'log');
