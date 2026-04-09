function [Ex_out,Ey_out,Ez_out,Hx_out,Hy_out,Hz_out] = scattering(r_eff,gird,Lambda,Np_shape,Structure,IB,E01,H01,K01,CT,GPU,np,nb,plane_2D)
%处理相对介电常数
eps_nps=(np).^2;%粒子的
epsb=nb^2;%环境
if GPU==1
    eps_NP_eb=gpuArray(eps_nps./epsb);  % Ratio of metal-to-medium dielectric function
elseif GPU==0
    eps_NP_eb=eps_nps./epsb;  % Ratio of metal-to-medium dielectric function
end
%直径
d_eff=2*r_eff;
%最小分割单元的尺寸
d=r_eff/gird;
%GPU适配
A0 = 1e5;
if GPU==1 
    E0=gpuArray(E01)*A0;
    H0=gpuArray(H01)*A0;
    K0=gpuArray(K01);   
elseif GPU==0
    E0=E01*A0;               % Incident electric field
    H0=H01*A0;
    K0=K01;               % unit vector in direction of wave vector
end
%波矢
k=2*pi./Lambda*nb;
kvec=k*K0;
%粒子结构处理
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
%入射光场处理
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
%找到包络粒子的矩形
[Lx,Ly,Lz] = Nps_parameters(r_eff,Np_shape);
%xyz方向最大的长度，各个方向上偶极子的数量和总的偶极子的数量
[Max_x,Max_y,Max_z,N,Nx,Ny,Nz,r_block,X,Y,Z,d_inter]=Coordinates(GPU,d,Lx,Ly,Lz,...
    d_eff,Structure,arrangement); 
%给出偶极子的位置索引
[INDEX_INSIDE]=INDEX_INSIDE_NP(GPU,X,Y,Z,N,Np_shape,Lx,Ly,Lz,Structure,d_inter,arrangement);
INDEX_IN=reshape(INDEX_INSIDE,[Nx,Ny,Nz]);
%计算相互作用矩阵函数的组成部分
[rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,rjkrjk31_I,rjkrjk32_I,...
    rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,RJK]=RijRij(r_block);
%计算极化率（的倒数/逆）
[Inverse_Alpha]=Polarizability(GPU,kvec,eps_NP_eb,INDEX_IN,d,E0);
%给出入射电场
[E_x,E_y,E_z,E_vector]=Incident_Field(Lambda,IB,nb,r_block,kvec,K0,INDEX_INSIDE,Nx,Ny,Nz,E0,z0,Waist_r,p,l);
%相互作用矩阵表达式中的系数
Exp_ikvec_rjk=exp(1i*norm(kvec)*RJK)./RJK;
ikvec_rjk=(1i*norm(kvec)*RJK-1)./(RJK.^2);
%计算相互作用矩阵
[Axx,Axy,Axz,Ayy,Ayz,Azz]=Interaction_Matrix(kvec,Exp_ikvec_rjk,...
    ikvec_rjk, rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,...
    rjkrjk31_I,rjkrjk32_I,rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,Nx,Ny,Nz);
%傅里叶变换
[FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ]=FFT_Interaction(GPU,Axx...
    ,Axy,Axz,Ayy,Ayz,Azz,Nx,Ny,Nz);
%共轭梯度下降，得到偶极子强度
[px,py,pz]=Biconjugate_Gradient(E_x,E_y,E_z,Nx,Ny,Nz,N,Inverse_Alpha,...
    INDEX_IN,E_vector,FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ,CT);
px=px.*INDEX_IN;
py=py.*INDEX_IN;
pz=pz.*INDEX_IN;
%散射平面的参数
[r_block,r_plane,Nx,Ny,Nz,N_plane,x_plane,y_plane,z_plane,SIZE,Xex,Yex,Zex]=...
    Two_D_plane_for_forwardscattering(GPU,d,Max_x,Max_y,Max_z,Lx,Ly,Lz,plane_2D);
%入射平面上的电磁场
[Ex_incident,Ey_incident,Ez_incident,Hx_incident,Hy_incident,Hz_incident]=Incident_Field_plane(r_plane,E0,H0,K0,z0,kvec,Waist_r,IB,Lambda,nb,p,l);
%粒子坐标
r_block_forEH=block(GPU,d,Max_x,Max_y,Max_z);
%计算电磁场
[Ex_out, Ey_out, Ez_out, Hx_out, Hy_out, Hz_out] = EM_out(r_block_forEH, r_plane, N_plane, px, py, pz, kvec, Ex_incident, Ey_incident, Ez_incident, Hx_incident, Hy_incident, Hz_incident);
%单独处理穿过粒子的截面
[Outside_Index]=Excluding_NPs(x_plane,y_plane,z_plane,N_plane,Np_shape,Lx,Ly,Lz,d_inter,Structure,arrangement);
if plane_2D=="yz"||plane_2D=="xz"||plane_2D=="xy"
    Ex_out = Ex_out.*Outside_Index;
    Ey_out = Ey_out.*Outside_Index;
    Ez_out = Ez_out.*Outside_Index;
    Hx_out = Hx_out.*Outside_Index;
    Hy_out = Hy_out.*Outside_Index;
    Hz_out = Hz_out.*Outside_Index;
end
%整理格式
Ex_out=reshape(Ex_out,[SIZE(1,1),SIZE(1,2)]);
Ey_out=reshape(Ey_out,[SIZE(1,1),SIZE(1,2)]);
Ez_out=reshape(Ez_out,[SIZE(1,1),SIZE(1,2)]);
Hx_out=reshape(Hx_out,[SIZE(1,1),SIZE(1,2)]);
Hy_out=reshape(Hy_out,[SIZE(1,1),SIZE(1,2)]);
Hz_out=reshape(Hz_out,[SIZE(1,1),SIZE(1,2)]);

