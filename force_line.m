function [F_line,T_line] = force_line(r_eff,gird,Lambda,Np_shape,Structure,IB,E01,H01,K01,CT,GPU,np,nb,direction,l_scan,z0,p,l,Waist_r)
eps_nps=(np).^2;%粒子的
epsb=nb^2;%环境
if GPU==1
    eps_NP_eb=gpuArray(eps_nps./epsb);  % Ratio of metal-to-medium dielectric function
elseif GPU==0
    eps_NP_eb=eps_nps./epsb;  % Ratio of metal-to-medium dielectric function
end
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

%找到包络粒子的矩形
[Lx, Ly, Lz] = Nps_parameters(r_eff, Np_shape);
%xyz方向最大的长度，各个方向上偶极子的数量和总的偶极子的数量
[Max_x, Max_y, Max_z, N, Nx, Ny, Nz, r_block_orig, X, Y, Z, d_inter] = Coordinates(...
    GPU, d, Lx, Ly, Lz, 2*r_eff, Structure, arrangement);
%给出偶极子的位置索引
[INDEX_INSIDE] = INDEX_INSIDE_NP(GPU, X, Y, Z, N, Np_shape, Lx, Ly, Lz, Structure, d_inter, arrangement);
INDEX_IN = reshape(INDEX_INSIDE, [Nx, Ny, Nz]);     % 三维索引矩阵
%计算相互作用矩阵函数的组成部分
[rjkrjk1_I, rjkrjk2_I, rjkrjk3_I, rjkrjk4_I, rjkrjk5_I, rjkrjk6_I, ...
    rjkrjk31_I, rjkrjk32_I, rjkrjk33_I, rjkrjk34_I, rjkrjk35_I, rjkrjk36_I, RJK] = RijRij(r_block_orig);
%遍历位置和循环预处理
N_scan = length(l_scan);                           % 遍历点数
if GPU == 1
    F_line = zeros(N_scan, 3, 'gpuArray');      % 梯度力存储（N_scan×3）
    T_line = zeros(N_scan, 3, 'gpuArray');     % 力矩存储
else
    F_line = zeros(N_scan, 3);
    T_line = zeros(N_scan, 3);
end
%遍历位置信息
for i = 1:N_scan
    fprintf('位置 %d\n', l_scan(i));
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
    %入射电场
    [E_x, E_y, E_z, E_vector] = Incident_Field(...
        Lambda, IB, nb, r_block, kvec, K0, INDEX_INSIDE, Nx, Ny, Nz, E0, z0, Waist_r,p,l);
    %计算极化率（的倒数/逆）
    [Inverse_Alpha] = Polarizability(GPU, kvec, eps_NP_eb, INDEX_IN, d, E0);
    %相互作用矩阵表达式中的系数
    Exp_ikvec_rjk = exp(1i * norm(kvec) * RJK) ./ RJK;
    ikvec_rjk = (1i * norm(kvec) * RJK - 1) ./ (RJK.^2);
    %计算相互作用矩阵
    [Axx, Axy, Axz, Ayy, Ayz, Azz] = Interaction_Matrix(kvec, Exp_ikvec_rjk, ikvec_rjk, rjkrjk1_I, rjkrjk2_I, rjkrjk3_I,rjkrjk4_I, rjkrjk5_I, rjkrjk6_I, rjkrjk31_I, rjkrjk32_I, rjkrjk33_I,rjkrjk34_I, rjkrjk35_I, rjkrjk36_I, Nx, Ny, Nz);
    %傅里叶变换
    [FFT_AXX, FFT_AXY, FFT_AXZ, FFT_AYY, FFT_AYZ, FFT_AZZ] = FFT_Interaction(GPU, Axx, Axy, Axz, Ayy, Ayz, Azz, Nx, Ny, Nz);
    %共轭梯度下降，得到偶极子强度
    [px, py, pz] = Biconjugate_Gradient(E_x, E_y, E_z, Nx, Ny, Nz, N, Inverse_Alpha, INDEX_IN, E_vector,FFT_AXX, FFT_AXY, FFT_AXZ, FFT_AYY, FFT_AYZ, FFT_AZZ, CT);
    px = px .* INDEX_IN;
    py = py .* INDEX_IN;
    pz = pz .* INDEX_IN;
    % [F_line(i, :),T_line(i, :)] = scattering_force_P2F(px,py,pz,arrangement,GPU,d,Max_x,Max_y,Max_z,Lx,Ly,Lz,x_pos,y_pos,z_pos, ...
    %     IB,K01,Waist_r,Lambda,nb,z0,E01,H01,INDEX_IN,Inverse_Alpha,Nx,Ny,Nz,Np_shape,kvec,d_inter,Structure,N,p,l);
    [F_line(i, :),T_line(i, :)] = scattering_force_P2F(px,py,pz,GPU,d,Max_x,Max_y,Max_z,Lx,Ly,Lz,x_pos,y_pos,z_pos,IB,K0,Waist_r,Lambda,nb,z0,E0,H0,INDEX_IN,kvec,p,l);
end
end

