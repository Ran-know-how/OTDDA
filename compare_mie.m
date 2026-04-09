clear all
close all
clc
% %% DDA
parallel.gpu.enableCUDAForwardCompatibility(true)
CT=10^(-10);% Convergence threshold value in interative solver迭代求解器中的收敛阈值
%============ LSPR wavelength and bulk RI of metal at LSPR ===============%
%=========================================================================%
GPU=1;
clc
Lambda=532e-9;
Re_n=1.5;
Im_n=0;
clc
eps=(Re_n+1i*Im_n).^2;
Re_eps=real(eps);
Im_eps=imag(eps);
nb=1;
epsb=nb^2;         % Dielectric function of background medium
eps0 = 8.854187817e-12; % 真空介电常数
mu0 = 4*pi*1e-7;   
r_eff=1000e-9;
d_eff=2*r_eff;
volume=4*pi/3*(r_eff^3);
d=r_eff/20;
% Np_shape='ellipsoid'; 
% Np_shape= 'rod';
% Np_shape= 'rec_block';
Np_shape= 'spherical';
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
[Lx,Ly,Lz] = Nps_parameters(r_eff,Np_shape);%返回包络住粒子的矩形，用矩形是因为FFT需要用矩形
[Max_x,Max_y,Max_z,N,Nx,Ny,Nz,r_block,X,Y,Z,d_inter]=Coordinates(GPU,d,Lx,Ly,Lz,...
    d_eff,Structure,arrangement); %返回了xyz方向最大的长度，各个方向上偶极子的数量和总的偶极子的数量
Nx_target=Nx;
Ny_target=Ny; 
Nz_target=Nz;
Xrectan=X;
Yrectan=Y;
Zrectan=Z;
c0 = 3e8;
E01=[1 0 0];%偏振方向为x
H01=[0 1 0]/(120*pi);
K01=[0 0 1];%传播方向为z正方向
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
% IB="gaussian" ; 
IB="plane wave" ;
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
    l= 0;
elseif IB=="laguerre-gaussian"
    z0=0;              % Focus point of the Gaussian beam
    NA = 0.3;
    n_object = 1;
    Waist_r = 1/(pi*NA*sqrt(n_object));     
    p = 0;
    l= 0;
end 
[INDEX_INSIDE]=INDEX_INSIDE_NP(GPU,X,Y,Z,N,Np_shape,Lx,Ly,Lz,Structure,d_inter,arrangement);
%INDEX_INSIDE格式是[m,1]，m是偶极子数量，（x，1）中存放的值是第m个坐标位置是偶极子
INDEX_IN=reshape(INDEX_INSIDE,[Nx,Ny,Nz]);% INDEX_IN contains zeros and ones, 
                  % ... zeros for dipole outside of NP, and ones for inside    转化为一列便于计算                                 
%=========================================================================%
[rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,rjkrjk31_I,rjkrjk32_I,...
    rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,RJK]=RijRij(r_block);
eps_NP_eb=ep_nps_eb;
kvec=k*K0;
[Inverse_Alpha]=Polarizability(GPU,kvec,eps_NP_eb,INDEX_IN,d,E0);
[E_x,E_y,E_z,E_vector]=Incident_Field(Lambda,IB,nb,r_block,kvec,K0,INDEX_INSIDE,Nx,Ny,Nz,E0,z0,Waist_r,p,l);
Exp_ikvec_rjk=exp(1i*norm(kvec)*RJK)./RJK;%A表达式中的一个系数
ikvec_rjk=(1i*norm(kvec)*RJK-1)./(RJK.^2);%A表达式中的一个系数   % ikvec_rjk=(1i*norm(kvec)*rjk-1)/rjk^2
[Axx,Axy,Axz,Ayy,Ayz,Azz]=Interaction_Matrix(kvec,Exp_ikvec_rjk,...
    ikvec_rjk, rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,...
    rjkrjk31_I,rjkrjk32_I,rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,Nx,Ny,Nz);
[FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ]=FFT_Interaction(GPU,Axx...
    ,Axy,Axz,Ayy,Ayz,Azz,Nx,Ny,Nz);
%计算了三个方向上相互作用矩阵的傅里叶变化的值，需要注意的是对其进行了偶扩展，三个方向上面的扩展顺序
[px,py,pz]=Biconjugate_Gradient(E_x,E_y,E_z,Nx,Ny,Nz,N,Inverse_Alpha,...
    INDEX_IN,E_vector,FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ,CT);
%=========================================================================%
clear FFT_AXX FFT_AXY FFT_AXZ FFT_AYY FFT_AYZ FFT_AZZ 
%=========================================================================%
%非偶极子索引的地方被舍去了，值得研究      
px=px.*INDEX_IN;
py=py.*INDEX_IN;
pz=pz.*INDEX_IN;
PX_vector=reshape(px,[N,1]);
PY_vector=reshape(py,[N,1]);
PZ_vector=reshape(pz,[N,1]); 
clear  a_CM_Nps anr_Nps aLDR_Nps  a_CM_Matrix anr_Matrix aLDR_Matrix ...
    Ex Ey Ez E_x E_y E_z ...
    E_vector Exp_ikvec_rjk ikvec_rjk r_block 
plane_2D="yz_x=2N";
% plane_2D="yz_x=-2N";
% plane_2D="xz_y=2N";
% plane_2D="xz_y=-2N";
% plane_2D="xy_z=2N";
% plane_2D="xy_z=-2N";
% plane_2D="xy" ;
% plane_2D="xz";
% plane_2D="yz"; 
[r_block,r_plane,Nx,Ny,Nz,N_plane,x_plane,y_plane,z_plane,SIZE,Xex,Yex,Zex]=...
    Two_D_plane_for_forwardscattering(GPU,d,Max_x,Max_y,Max_z,Lx,Ly,Lz,plane_2D);
[Ex_incident,Ey_incident,Ez_incident,Hx_incident,Hy_incident,Hz_incident]=Incident_Field_plane(r_plane,E0,H0,K0,z0,kvec,Waist_r,IB,Lambda,nb,p,l);
[r_block_forH]=block(GPU,d,Max_x,Max_y,Max_z);
[Ex_out, Ey_out, Ez_out, Hx_out, Hy_out, Hz_out] = EM_out(r_block_forH, r_plane, N_plane, px, py, pz, kvec, Ex_incident, Ey_incident, Ez_incident, Hx_incident, Hy_incident, Hz_incident);
Ex_out=reshape(Ex_out,[SIZE(1,1),SIZE(1,2)]);
Ey_out=reshape(Ey_out,[SIZE(1,1),SIZE(1,2)]);
Ez_out=reshape(Ez_out,[SIZE(1,1),SIZE(1,2)]);
Hx_out=reshape(Hx_out,[SIZE(1,1),SIZE(1,2)]);
Hy_out=reshape(Hy_out,[SIZE(1,1),SIZE(1,2)]);
Hz_out=reshape(Hz_out,[SIZE(1,1),SIZE(1,2)]);
Sx = real((Ey_out.*conj(Hz_out)-Ez_out.*conj(Hy_out)));
Sy = real((Ez_out.*conj(Hx_out)-Ex_out.*conj(Hz_out)));
Sz = real((Ex_out.*conj(Hy_out)-Ey_out.*conj(Hx_out)));
E_DDA = sqrt(Ex_out.*conj(Ex_out)+Ey_out.*conj(Ey_out)+Ez_out.*conj(Ez_out));
H_DDA = sqrt(Hx_out.*conj(Hx_out)+Hy_out.*conj(Hy_out)+Hz_out.*conj(Hz_out));
S_DDA = sqrt(Sx.*conj(Sx)+Sy.*conj(Sy)+Sz.*conj(Sz));
figure;
subplot(331)
imagesc(abs(Ex_out))
title('Ex')
axis off
subplot(332)
imagesc(abs(Ey_out))
title('Ey')
axis off
subplot(333)
imagesc(abs(Ez_out))
title('Ez')
axis off
subplot(334)
imagesc(abs(Hx_out))
title('Hx')
axis off
subplot(335)
imagesc(abs(Hy_out))
title('Hy')
axis off
subplot(336)
imagesc(abs(Hz_out))
title('Hz')
axis off
subplot(337)
imagesc(Sx)
title('Sx')
axis off
subplot(338)
imagesc(Sy)
title('Sy')
axis off
subplot(339)
imagesc(Sz)
title('Sz')
axis off
%% mie
%沿z方向传播
m1 = 1.5;%(*颗粒球体折射率（二氧化硅1.45-1.5，润滑油1.4左右）*)
m2 = 1  ;%(*周围介质的折射率（空气1）*)
M = m1/m2;%(*颗粒球体对周围介质的相对折射率*)
a = r_eff;
Lambda0 = 532*10^-9;%(*激光在真空中的波长*)
Lambda = Lambda0/m2;%(*激光在介质中的波长*)
k = 2*pi*Lambda;%(*激光在介质中的波数*)
x = k*a;%(*尺寸参数*)
Epsilon0 = 8.85*10^-12;%(*真空介电常数*)
Mu0 = 4*pi*10^-7;%(*真空磁导率*)
ns = round(x + 4*(x^(1/3)) + 2);%(*n的最大取值*)
screen_resolution = 201; % 屏幕分辨率
% screen_distance =0; % 屏幕距离
% screen_distance = 1e-8; % 屏幕距离
screen_distance =a*2; % 屏幕距离
screen_size = 2*a;% 屏幕大小
% % %前向散射
% [xc, yc] = meshgrid(linspace(-screen_size, screen_size, screen_resolution), linspace(-screen_size, screen_size, screen_resolution));
% zc = ones(size(xc)) * screen_distance;
% 侧向散射
[yc, zc] = meshgrid(linspace(-screen_size, screen_size, screen_resolution), linspace(-screen_size, screen_size, screen_resolution));
xc = ones(size(zc)) * screen_distance;
% [xc, zc] = meshgrid(linspace(-screen_size, screen_size, screen_resolution), linspace(-screen_size, screen_size, screen_resolution));
% yc = ones(size(zc)) * screen_distance;
[E,H] = calcmie_nf(a,m1,m2,Lambda,xc,yc,zc);
%这里为了和DDA对应进行了坐标变换
Ex = E(:,:,1);
Ey = E(:,:,2);
Ez = E(:,:,3);
Hx = H(:,:,1);
Hy = H(:,:,2);
Hz = H(:,:,3);
%加上入射场
Sx = 0.5*real(Ey.*conj(Hz) - Ez.*conj(Hy));
Sy = 0.5*real(Ez.*conj(Hx) - Ex.*conj(Hz));
Sz = 0.5*real(Ex.*conj(Hy) - Ey.*conj(Hx));
E_mie = sqrt(Ex.*conj(Ex)+Ey.*conj(Ey)+Ez.*conj(Ez));
H_mie = sqrt(Hy.*conj(Hy)+Hx.*conj(Hx)+Hz.*conj(Hz));
S_mie = sqrt(Sx.*conj(Sx)+Sy.*conj(Sy)+Sz.*conj(Sz));
figure;
title("Mie")
subplot(331)
imagesc(abs(Ex))
title('Ex')
axis off
subplot(332)
imagesc(abs(Ey))
title('Ey')
axis off
subplot(333)
imagesc(abs(Ez))
title('Ez')
axis off
subplot(334)
imagesc(abs(Hx))
title('Hx')
axis off
subplot(335)
imagesc(abs(Hy))
title('Hy')
axis off
subplot(336)
imagesc(abs(Hz))
title('Hz')
axis off
subplot(337)
imagesc(Sx)
title('Sx')
axis off
subplot(338)
imagesc(Sy)
title('Sy')
axis off
subplot(339)
imagesc(Sz)
title('Sz')
axis off
% figure;
% subplot(231)
% imagesc(abs(E_DDA))
% subplot(232)
% imagesc(abs(H_DDA))
% subplot(233)
% imagesc(abs(S_DDA))
% subplot(234)
% imagesc(abs(E_mie))
% subplot(235)
% imagesc(abs(H_mie))
% subplot(236)
% imagesc(abs(S_mie))
% figure('Position', [100, 100, 200, 200]);
% imagesc(abs(E_DDA))
% axis off
% figure('Position', [100, 100, 200, 200]);
% imagesc(abs(E_mie))
% axis off
