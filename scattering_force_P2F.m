function [F_total,N_total] = scattering_force_P2F(px,py,pz,GPU,d,Max_x,Max_y,Max_z,Lx,Ly,Lz,x_pos,y_pos,z_pos,IB,K0,Waist_r,Lambda,nb,z0,E0,H0,INDEX_IN,kvec,p,l)
eps0 = 8.854187817e-12; % 真空介电常数
mu0 = 4*pi*1e-7;        % 真空磁导率
px=px.*INDEX_IN;
py=py.*INDEX_IN;
pz=pz.*INDEX_IN;

plane_list = ["yz_x=-2N", "yz_x=2N", "xz_y=-2N", "xz_y=2N", "xy_z=-2N", "xy_z=2N"];
num_planes = length(plane_list);

plane_data = struct();
for i = 1:num_planes
    plane_2D = plane_list(i);

[r_block,r_plane,Nx,Ny,Nz,N_plane,x_plane,y_plane,z_plane,SIZE,Xex,Yex,Zex]=...
    Two_D_plane_for_forwardscattering(GPU,d,Max_x,Max_y,Max_z,Lx,Ly,Lz,plane_2D);
r_plane_inc = r_plane;

r_plane_inc(:,2) = r_plane_inc(:,2) + y_pos;   % y 偏移
r_plane_inc(:,3) = r_plane_inc(:,3) + z_pos;   % z 偏移
r_plane_inc(:,1) = r_plane_inc(:, 1) + x_pos;  % x分量偏移，实现粒子移动

[Ex_incident,Ey_incident,Ez_incident,Hx_incident,Hy_incident,Hz_incident]=Incident_Field_plane(r_plane_inc,E0,H0,K0,z0,kvec,Waist_r,IB,Lambda,nb,p,l);
[r_block_forH]=block(GPU,d,Max_x,Max_y,Max_z);
[Ex_out, Ey_out, Ez_out, Hx_out, Hy_out, Hz_out] = EM_out(r_block_forH, r_plane, N_plane, px, py, pz, kvec, Ex_incident, Ey_incident, Ez_incident, Hx_incident, Hy_incident, Hz_incident);
    plane_data(i).Ex = Ex_out;
    plane_data(i).Ey = Ey_out;
    plane_data(i).Ez = Ez_out;
    plane_data(i).Hx = Hx_out;
    plane_data(i).Hy = Hy_out;
    plane_data(i).Hz = Hz_out;
    plane_data(i).X = Xex;
    plane_data(i).Y = Yex;
    plane_data(i).Z = Zex;
    plane_data(i).normal = get_normal_vector(plane_2D); % 获取法向量
    plane_data(i).dA = d * d * ones(size(Xex)); % 计算面积元
end
F_total = zeros(1, 3);
N_total = zeros(1, 3);
for i = 1:num_planes
    
    %% ====== 1. 提取本平面的场和几何信息 ======
    Ex = plane_data(i).Ex(:);
    Ey = plane_data(i).Ey(:);
    Ez = plane_data(i).Ez(:);
    Hx = plane_data(i).Hx(:);
    Hy = plane_data(i).Hy(:);
    Hz = plane_data(i).Hz(:);
    % 复数场允许存在
    normal = plane_data(i).normal(:).';   % [nx,ny,nz]
    dA = plane_data(i).dA(:);             % 面积元

    X = plane_data(i).X(:);
    Y = plane_data(i).Y(:);
    Z = plane_data(i).Z(:);

    %% ====== 2. 计算时均应力张量 T_ij = (1/2) Re( eps0 Ei Ej* + mu0 Hi Hj* - 1/2(|E|^2 + |H|^2)δij ) ======
    absE2 = abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;
    absH2 = abs(Hx).^2 + abs(Hy).^2 + abs(Hz).^2;

    % 对角分量
    T_xx = 0.5 * real( eps0*(Ex .* conj(Ex)) + mu0*(Hx .* conj(Hx)) ...
        - 0.5*(eps0*absE2 + mu0*absH2) );
    T_yy = 0.5 * real( eps0*(Ey .* conj(Ey)) + mu0*(Hy .* conj(Hy)) ...
        - 0.5*(eps0*absE2 + mu0*absH2) );
    T_zz = 0.5 * real( eps0*(Ez .* conj(Ez)) + mu0*(Hz .* conj(Hz)) ...
        - 0.5*(eps0*absE2 + mu0*absH2) );

    % 非对角分量 T_ij = T_ji
    T_xy = 0.5 * real( eps0*(Ex .* conj(Ey)) + mu0*(Hx .* conj(Hy)) );
    T_xz = 0.5 * real( eps0*(Ex .* conj(Ez)) + mu0*(Hx .* conj(Hz)) );
    T_yz = 0.5 * real( eps0*(Ey .* conj(Ez)) + mu0*(Hy .* conj(Hz)) );

    %% ====== 3. traction: F_density = T · n  ======
    nx = normal(1); ny = normal(2); nz = normal(3);

    F_density_x = T_xx*nx + T_xy*ny + T_xz*nz;
    F_density_y = T_xy*nx + T_yy*ny + T_yz*nz;
    F_density_z = T_xz*nx + T_yz*ny + T_zz*nz;

    %% ====== 4. 力元：dF = traction * dA ======
    dFx = F_density_x .* dA;
    dFy = F_density_y .* dA;
    dFz = F_density_z .* dA;

    %% ====== 5. 本平面的总力 ======
    F_plane = [sum(dFx), sum(dFy), sum(dFz)];
    %% ====== 6. 计算质心（每个平面独立） ======
    r_center = [mean(X), mean(Y), mean(Z)];

    %% ====== 7. 力矩：τ = r × dF ======
    Rx = X - r_center(1);
    Ry = Y - r_center(2);
    Rz = Z - r_center(3);

    Tx = Ry .* dFz - Rz .* dFy;
    Ty = Rz .* dFx - Rx .* dFz;
    Tz = Rx .* dFy - Ry .* dFx;

    N_plane = [sum(Tx), sum(Ty), sum(Tz)];

    %% ====== 8. 累计到 total ======
    F_total = F_total + F_plane;
    N_total = N_total + N_plane;
end
F_total = -F_total;
N_total = -N_total;%平移的是入射电场坐标，换成移动粒子应该是相对负值
end
