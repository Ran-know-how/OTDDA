%计算了粒子内部的电场，并以此为局部热源计算了颗粒内部的温度分布，画出了表面的温度分布情况
%代码应该正确
clear 
close all
clc
parallel.gpu.enableCUDAForwardCompatibility(true)
CT=10^(-5);                     % Convergence threshold value in interative solver迭代求解器中的收敛阈值

%============ LSPR wavelength and bulk RI of metal at LSPR ===============%
%=========================================================================%
GPU=1;
clc
Lambda=532e-9;
% Re_n=1.37;
% Im_n=7.57;
Re_n = 1.5;
Im_n = 0.1;
n = Re_n+Im_n;
eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;
clc

eps=(Re_n+1i*Im_n).^2;
Re_eps=real(eps);
Im_eps=imag(eps);

nb=1;
epsb=nb^2;         % Dielectric function of background medium
clc

r_eff=500e-9;
d_eff=2*r_eff;
volume=4*pi/3*(r_eff^3);
clc
d=r_eff/10;
clc
% Np_shape='ellipsoid'; 
Np_shape= 'rod';
% Np_shape= 'rec_block';
% Np_shape= 'spherical';
clc
Structure="monomeric";
% Structure="dimeric";
if Structure=="dimeric"
    fprintf('\nThree kind of arragments for neighboring NPs:');
    fprintf('\n1. For head to tail orientation in z-direction, type "z_orient"  \n');
    fprintf('2. For side by side orientation in x-direction, type "x_orient"  \n');
    fprintf('3. For side by side orientation in y-direction,type "y_orient"  \n');
    arrangement=input('\nType the orientation of the neighboring NPs:');
    clc
else
     arrangement="monomeric";
end
clc
r_np=r_eff;               % effective radius of the nanoparticle纳米颗粒的有效半径
W=2*pi*3*10^8./Lambda;         % incident radiation frequency入射辐射频率
omega = W;
k=2*pi./Lambda*nb;              % wave vector of light in first layer第一层中光的波矢量
eps_nps=eps; % Modified dielectric function
Refractive=eps_nps.^(0.5);
Re_Refractive=real(Refractive);
Im_Refractive=imag(Refractive); 

if GPU==1
    ep_nps_eb=gpuArray(eps_nps./epsb);  % Ratio of metal-to-medium dielectric function
elseif GPU==0
    ep_nps_eb=eps_nps./epsb;  % Ratio of metal-to-medium dielectric function
end
%=========================================================================%
%====== Finding dimension of the Rectangular block emcopasses NPs ========%
%=========================================================================%
[Lx,Ly,Lz] = Nps_parameters(r_eff,Np_shape);%返回包络住粒子的矩形，用矩形是因为FFT需要用矩形
clc
%=========================================================================%
%==== Finding coordinate of the dipoles inside the rectangular block =====%
%=========================================================================%
[Max_x,Max_y,Max_z,N,Nx,Ny,Nz,r_block,X,Y,Z,d_inter]=Coordinates(GPU,d,Lx,Ly,Lz,...
    d_eff,Structure,arrangement); %返回了xyz方向最大的长度，各个方向上偶极子的数量和总的偶极子的数量
Nx_target=Nx;
Ny_target=Ny;
Nz_target=Nz;
Xrectan=X;
Yrectan=Y;
Zrectan=Z;
%=========================================================================%

%===================== Selecting the incident light ======================%
%=========================================================================%
c0 = 3e8;
E0 = 1e5;
E01=[0 0 1]*E0;%偏振方向为z
H01=[0 1 0]*E0/(120*pi);
K01=[1 0 0];%传播方向为x
%==============%
if GPU==1 
    E0=gpuArray(E01);
    H0=gpuArray(H01);
    K0=gpuArray(K01);   
elseif GPU==0
    E0=E01;               % Incident electric field
    H0=H01;
    K0=K01;               % unit vector in direction of wave vector
end
%=============================%
clc
IB="gaussian" ; 
% IB="plane wave" ;
% IB="laguerre-gaussian" ;
clc
if IB=="plane wave" 
    z0=0;              % Focus point of the Gaussian beam
    Waist_r=100;       % ratio of waist raduis of beam to wavelength
    p = 0;
    l= 0;
elseif IB=="gaussian"
    z0=0;              % Focus point of the Gaussian beam
    NA = 0.3;
    n_object = 1;
    Waist_r = 1/(pi*NA*sqrt(n_object));  
    p = 0;
    l = 0;
elseif IB=="laguerre-gaussian"
    z0=0;              % Focus point of the Gaussian beam
    NA = 0.7;
    n_object = 1;
    Waist_r = 1/(pi*NA*sqrt(n_object));     
    p = 2;
    l= 0;
end 

%=========================================================================%
% [INDEX_INSIDE]=INDEX_INSIDE_NP(GPU,X,Y,Z,N,Np_shape,Lx,Ly,Lz,Structure,d_inter,arrangement,theta);
[INDEX_INSIDE]=INDEX_INSIDE_NP(GPU,X,Y,Z,N,Np_shape,Lx,Ly,Lz,Structure,d_inter,arrangement);
%INDEX_INSIDE格式是[m,1]，m是偶极子数量，（x，1）中存放的值是第m个坐标位置是偶极子
INDEX_IN=reshape(INDEX_INSIDE,[Nx,Ny,Nz]);% INDEX_IN contains zeros and ones, 
                  % ... zeros for dipole outside of NP, and ones for inside    转化为一列便于计算                                 
%=========================================================================%
[rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,rjkrjk31_I,rjkrjk32_I,...
    rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,RJK]=RijRij(r_block);
%rjkrjk中有6个组成部分（矩阵），其余位置为0
%=========================================================================%
eps_NP_eb=ep_nps_eb;
kvec=k*K0;
[Inverse_Alpha]=Polarizability(GPU,kvec,eps_NP_eb,INDEX_IN,d,E0);
%=========================================================================%
%=========================================================================%
% r_block(:,2) = r_block(:,2) + 3e-7;   % z 偏移
[E_x,E_y,E_z,E_vector]=Incident_Field(Lambda,IB,nb,r_block,kvec,K0,INDEX_INSIDE,Nx,Ny,Nz,E0,z0,Waist_r,p,l);
%返回了每个点位上入射电场的值，E_vector为三个方向电场排列在一起的一个向量
%=========================================================================%
Exp_ikvec_rjk=exp(1i*norm(kvec)*RJK)./RJK;%A表达式中的一个系数
ikvec_rjk=(1i*norm(kvec)*RJK-1)./(RJK.^2);%A表达式中的一个系数   % ikvec_rjk=(1i*norm(kvec)*rjk-1)/rjk^2
[Axx,Axy,Axz,Ayy,Ayz,Azz]=Interaction_Matrix(kvec,Exp_ikvec_rjk,...
    ikvec_rjk, rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,...
    rjkrjk31_I,rjkrjk32_I,rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,Nx,Ny,Nz);
%把A的几个部分乘在一起算出来了，A为一个3*3的上三角矩阵
%=========================================================================%
[FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ]=FFT_Interaction(GPU,Axx...
    ,Axy,Axz,Ayy,Ayz,Azz,Nx,Ny,Nz);
%计算了三个方向上相互作用矩阵的傅里叶变化的值，需要注意的是对其进行了偶扩展，三个方向上面的扩展顺序
%=========================================================================%
%=========================================================================%
[px,py,pz] = Biconjugate_Gradient(E_x,E_y,E_z,Nx,Ny,Nz,N,Inverse_Alpha,...
    INDEX_IN,E_vector,FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ,CT);%原本的迭代方法
              
clear FFT_AXX FFT_AXY FFT_AXZ FFT_AYY FFT_AYZ FFT_AZZ 
%=========================================================================%
%粒子之外的极化强度设定为0
px=px.*INDEX_IN;
py=py.*INDEX_IN;
pz=pz.*INDEX_IN;

PX_vector=reshape(px,[N,1]);
PY_vector=reshape(py,[N,1]);
PZ_vector=reshape(pz,[N,1]);
P_vector=[PX_vector;PY_vector;PZ_vector];
%=========================================================================%   
                    % Deleting unnessesary data %
%=========================================================================%  
clear  a_CM_Nps anr_Nps aLDR_Nps  a_CM_Matrix anr_Matrix aLDR_Matrix ...
    Ex Ey Ez E_x E_y E_z ...
    E_vector Exp_ikvec_rjk ikvec_rjk 
%=========================================================================% 
Ex_whole = Inverse_Alpha.*px;
Ey_whole = Inverse_Alpha.*py;
Ez_whole = Inverse_Alpha.*pz;

E_whole = sqrt((Ex_whole).*conj(Ex_whole)+(Ey_whole).*conj(Ey_whole)+Ez_whole.*conj(Ez_whole));

epsilon_2 = imag(eps_nps);  % 这里的 eps_nps 已经定义为复介电函数

q_local = 0.5 * omega .* epsilon_2 .* (abs(E_whole).^2) * eps0;

a = r_eff;   % simply rename for consistency

% ---------------- thermal/gas inputs (与第二段完全一致) ----------------
k_p = 0.2;  
% k_p = 200;
k_g = 0.0257;         
lambda_mfp = 65e-9;   
beta = 0.8;           
T_inf = 300;          
Kn = lambda_mfp / a;
h_cont = k_g / a;
h_fm = k_g / (beta * lambda_mfp);
h_eff = (1 ./ (1./h_cont + 1./h_fm));
N = Nx;
L = a;
grid = linspace(-L, L, N);
dx = grid(2) - grid(1);
%% ======================================================
% ============ 1. 计算 3D 温度场 =========================
% =======================================================
disp("Solving for temperature field...")
T_sol = heat_solver_general( ...
        INDEX_IN, ...     % 颗粒内部
        q_local, ...      % 吸收功率密度（来自你的电场）
        dx, ...           % 网格步长
        k_p, ...          % 颗粒热导率
        h_eff, ...        % 有效对流（连续+自由分子）
        T_inf );          % 背景温度
%% ======================================================
% ============ 2. 计算表面 SDF + 法向 ====================
% =======================================================
disp("Computing surface normals...")
[surf_idx, nx, ny, nz] = surface_normals_3D(INDEX_IN, dx);

%% ======================================================
% ============ 3. 计算光泳力（Rohatschek 非球体）=========
% =======================================================
disp("Computing photophoretic force...")
[F_ph, M_ph] = photophoretic_force(a,surf_idx, nx, ny, nz, T_sol, X, Y, Z, dx, k_g, k_p, T_inf, Kn, beta);
fprintf("\nPhotophoretic Force F = [%.3e, %.3e, %.3e] N\n",F_ph(1),F_ph(2),F_ph(3));
fprintf("\nPhotophoretic Torque_total   = [%.3e, %.3e, %.3e] N\n",M_ph(1),M_ph(2),M_ph(3));
%% ======================================================
% ============ 4. 可视化表面温度 =========================
% =======================================================
% figure('Position', [100, 100, 800, 400]); hold on;
% subplot(121)
[sx,sy,sz] = ind2sub(size(T_sol),surf_idx);
figure('Position', [100, 100, 400, 350]);
scatter3(sx,sy,sz,20,T_sol(surf_idx),'filled');
colorbar; axis equal
xlabel("x"); ylabel("y"); zlabel("z");
title("Surface temperature distribution")
view(45,25);
% subplot(122)
figure('Position', [100, 100, 400, 350]);
scatter3(sx,sy,sz,20,abs(E_whole(surf_idx)),'filled');
colorbar; axis equal
xlabel("x"); ylabel("y"); zlabel("z");
title("Surface E distribution")
view(45,25);
