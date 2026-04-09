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
r_eff=1000e-9;
%分割
grid = 25;
%粒子形状：仅为球形
Np_shape= 'spherical';
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
[T_sol,E_whole,inside_idx,XX,YY,ZZ] = temperature(r_eff,grid...
    ,Lambda,Np_shape,Structure,IB,E01,H01,K01,CT,GPU,np,nb);

figure('Name','Temperature and |E|^2 inside sphere','Position',[100 100 1200 500]);
subplot(1,2,1);
scatter3(XX(inside_idx)*1e6, YY(inside_idx)*1e6, ZZ(inside_idx)*1e6, 20, T_sol(inside_idx), 'filled');
axis equal; colorbar; title('T inside (K)'); xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)');
view(45,25);

subplot(1,2,2);
scatter3(XX(inside_idx)*1e6, YY(inside_idx)*1e6, ZZ(inside_idx)*1e6, 20, E_whole(inside_idx), 'filled');
axis equal; colorbar; title('|E|^2 inside (V^2/m^2)'); xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)');
view(45,25);

%% ------------------ 输出 ----------------------
% fprintf('T_inside: min=%g, max=%g, mean=%g\n', ...
%     min(T_sol(mask)), max(T_sol(mask)), mean(T_sol(mask)));

[Fx,Fy,Fz] = Photophoretic_force(T_sol,XX,YY,ZZ,r_eff,grid);
F_vec = [Fx, Fy, Fz];

fprintf('Estimated photophoretic force (heuristic) = [%.3e, %.3e, %.3e] N\n', F_vec);
fprintf('Magnitude = %.3e N\n', norm(F_vec));                                                                                                                              