%遍历粒子的位置
%计算力的方法是T张量表面积分法

clc; clear all; close all;
parallel.gpu.enableCUDAForwardCompatibility(true); % GPU兼容性开启
CT = 10^(-5);                                     % 迭代收敛阈值
GPU = 1;                                           % GPU启用标识

% 2. 光镊核心参数（与T矩阵场景完全对齐）
% 2.1 激光与介质参数
Lambda = 532e-9;                                   % 真空激光波长（532nm，统一）
c0 = 3e8;
eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;
omega=2*pi*c0./Lambda;         % incident radiation frequency入射辐射频率
nb = 1;                                          % 背景折射率（空气，统一）                               % 真空光速（统一）
beam_power_ref = 1;                            % 参考激光功率（10mW，与T矩阵一致）

% 2.2 粒子参数
% Re_n=1.37;
% Im_n=7.57;
Re_n=1.5;
Im_n=0.1;
r_eff = 500e-9;                                    % 粒子有效半径（500nm，统一）
d = r_eff / 20;                                    % 偶极子间距（25nm，DDA精度足够）
Np_shape = 'spherical';                            % 球形粒子（统一）
Structure = "monomeric";                           % 单体结构（统一）
arrangement = "monomeric";                         % 单体排列

% 2.3 高斯光束参数（光镊核心，与T矩阵一致）
IB="gaussian" ; 
% IB="plane wave" ;
% IB="laguerre-gaussian" ;
clc
if IB=="plane wave" 
    z0_beam=0;              % Focus point of the Gaussian beam
    Waist_r=100;       % ratio of waist raduis of beam to wavelength
    p = 0;
    l= 0;
elseif IB=="gaussian"
    z0_beam=0;              % Focus point of the Gaussian beam
    NA = 0.3;
    n_object = 1;
    Waist_r = n_object/(pi*NA);  
    p = 0;
    l = 0;
elseif IB=="laguerre-gaussian"
    z0_beam=0;              % Focus point of the Gaussian beam
    NA = 0.5;
    n_object = 1;
    Waist_r = n_object/(pi*NA);     
    p = 4;
    l= 2;
end 
k = 2 * pi * nb / Lambda;                           % 介质中波矢

% 2.4 z方向位置遍历（与T矩阵完全相同：-2μm~+2μm，步长100nm）
direction = "y";
l_scan = -3*r_eff : r_eff/10 : 3*r_eff;              % 粒子中心z坐标遍历数组
N_scan = length(l_scan);                           % 遍历点数

% 3. 预计算DDA基础组件（仅计算1次，循环中复用）
% 3.1 粒子外包络矩形尺寸
[Lx, Ly, Lz] = Nps_parameters(r_eff, Np_shape);
% 3.2 偶极子坐标（原始，粒子中心在z=0）
[Max_x, Max_y, Max_z, N, Nx, Ny, Nz, r_block_orig, X, Y, Z, d_inter] = Coordinates(...
    GPU, d, Lx, Ly, Lz, 2*r_eff, Structure, arrangement);
% 3.3 粒子内部偶极子索引（粒子形状不变，索引固定）
[INDEX_INSIDE] = INDEX_INSIDE_NP(GPU, X, Y, Z, N, Np_shape, Lx, Ly, Lz, Structure, d_inter, arrangement);
INDEX_IN = reshape(INDEX_INSIDE, [Nx, Ny, Nz]);     % 三维索引矩阵
% 3.4 粒子介电参数（GPU适配）
eps = (Re_n + 1i*Im_n).^2;                         % 粒子介电函数
epsb = nb^2;                                       % 背景介电函数
eps_nps = eps;
if GPU == 1
    ep_nps_eb = gpuArray(eps_nps ./ epsb);         % 介电比（GPU）
else
    ep_nps_eb = eps_nps ./ epsb;                   % 介电比（CPU）
end
% 3.5 预计算偶极子间距离矩阵（RJK，仅1次）
[rjkrjk1_I, rjkrjk2_I, rjkrjk3_I, rjkrjk4_I, rjkrjk5_I, rjkrjk6_I, ...
    rjkrjk31_I, rjkrjk32_I, rjkrjk33_I, rjkrjk34_I, rjkrjk35_I, rjkrjk36_I, RJK] = RijRij(r_block_orig);

    E0 = 1e5;                                      % 电场振幅（可通过功率校准，此处与原代码一致）
    E01 = [0 0 1] * E0;                            % z方向偏振（与T矩阵偏振一致）
    K01 = [1 0 0];                                 % x方向传播（与T矩阵一致）
    if GPU == 1
        E0_gpu = gpuArray(E01);
        K0_gpu = gpuArray(K01);
        kvec = k * K0_gpu;
    else
        E0_gpu = E01;
        K0_gpu = K01;
        kvec = k * K0_gpu;
    end
% 4. 定义结果存储数组（GPU/CPU适配）
if GPU == 1
    F_vec = zeros(N_scan, 3, 'gpuArray');      % 梯度力存储（N_scan×3）
else
    F_vec = zeros(N_scan, 3);
end

% 5. 循环遍历z方向位置（核心：粒子在光镊中移动，计算每个位置受力）
for i = 1:N_scan
    % 5.1 当前粒子中心z位置（与T矩阵遍历位置一致）
    % z_pos = z_scan(i);
    r_block = r_block_orig;
    if direction=="z"
    z_pos = l_scan(i);
    y_pos = 0;
    x_pos = 0;
    r_block(:, 3) = r_block(:, 3) + z_pos;  % z分量偏移，实现粒子移动
    elseif direction=="y"
    y_pos = l_scan(i);
    z_pos = 0;
    x_pos = 0;
    r_block(:, 2) = r_block(:, 2) + y_pos;  % y分量偏移，实现粒子移动
    else 
    x_pos = l_scan(i);
    z_pos = 0;
    y_pos = 0;
    r_block(:, 1) = r_block(:, 1) + x_pos;  % y分量偏移，实现粒子移动
    end

    
    % 5.4 DDA核心计算（每个位置重新计算光场与偶极矩）
    % 5.4.1 入射电场（粒子在z_pos处的光场）
    [E_x, E_y, E_z, E_vector] = Incident_Field(...
        Lambda, IB, nb, r_block, kvec, K0_gpu, INDEX_INSIDE, Nx, Ny, Nz, E0_gpu, z0_beam, Waist_r,p,l);
    % 5.4.2 偶极子极化率
    [Inverse_Alpha] = Polarizability(GPU, kvec, ep_nps_eb, INDEX_IN, d, E0_gpu);
    % 5.4.3 相互作用矩阵与FFT
    Exp_ikvec_rjk = exp(1i * norm(kvec) * RJK) ./ RJK;
    ikvec_rjk = (1i * norm(kvec) * RJK - 1) ./ (RJK.^2);
    [Axx, Axy, Axz, Ayy, Ayz, Azz] = Interaction_Matrix(...
        kvec, Exp_ikvec_rjk, ikvec_rjk, rjkrjk1_I, rjkrjk2_I, rjkrjk3_I, ...
        rjkrjk4_I, rjkrjk5_I, rjkrjk6_I, rjkrjk31_I, rjkrjk32_I, rjkrjk33_I, ...
        rjkrjk34_I, rjkrjk35_I, rjkrjk36_I, Nx, Ny, Nz);
    [FFT_AXX, FFT_AXY, FFT_AXZ, FFT_AYY, FFT_AYZ, FFT_AZZ] = FFT_Interaction(...
        GPU, Axx, Axy, Axz, Ayy, Ayz, Azz, Nx, Ny, Nz);
    % 5.4.4 双共轭梯度求解偶极矩
    [px, py, pz] = Biconjugate_Gradient(...
        E_x, E_y, E_z, Nx, Ny, Nz, N, Inverse_Alpha, INDEX_IN, E_vector, ...
        FFT_AXX, FFT_AXY, FFT_AXZ, FFT_AYY, FFT_AYZ, FFT_AZZ, CT);
    
    % 5.5 筛选粒子内部偶极矩与总电场
    px = px .* INDEX_IN;
    py = py .* INDEX_IN;
    pz = pz .* INDEX_IN;
    Ex_whole = Inverse_Alpha .* px;
    Ey_whole = Inverse_Alpha .* py;
    Ez_whole = Inverse_Alpha .* pz;
    E_whole = sqrt((Ex_whole).*conj(Ex_whole)+(Ey_whole).*conj(Ey_whole)+Ez_whole.*conj(Ez_whole));
    clear E_x E_y E_z E_vector Inverse_Alpha Axx Axy Axz Ayy Ayz Azz FFT_*;
    F_vec(i,:) = heat_force_E2F(r_eff,E_whole,eps_nps,Lambda,Nx,Ny,Nz);
end

% 6. 结果格式转换（GPU→CPU，便于绘图）
if GPU == 1
    F_vec = gather(F_vec);
    l_scan = gather(l_scan);
end

figure('Position', [100, 100, 800, 400]);

plot(l_scan, F_vec);
legend Fx Fy Fz
xlabel('粒子y位置（m）');
ylabel('力（N）');
title('光镊纵向力（国际单位制）');
grid on;
% set(gca, 'YScale', 'log');


% E_z=abs(reshape(E_z,[SIZE(1,1),SIZE(1,2)]));
% imagesc(abs(E_z))
% title('Ez')
% axis off
