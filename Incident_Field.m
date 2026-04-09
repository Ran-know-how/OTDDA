function [E_x,E_y,E_z,E_vector]=Incident_Field(Lambda,IB,nb,r_block,kvec,K0,INDEX_INSIDE,Nx,Ny,Nz,E0,z0,Waist_r,p,l)
%输入变量：波长，入射波的类型：高斯波或者平面波，背景折射率，每个偶极子的坐标，波矢，波矢振动方向，偶极子索引，三个方向上的偶极子数量
%电场分量振动方向，高斯光斑的焦点位置，束腰半径与波长之比

kr=kvec(1)*r_block(:,1)+kvec(2)*r_block(:,2)+kvec(3)*r_block(:,3);
%kvec(1)kvec(2)kvec(3)分别是波矢的xyz三个方向分量
%波矢坐标分量乘以每个偶极子的xyz坐标
%得到了k·r的点乘，格式是[N,1]，kr（m,1）对应第m个偶极子的kr点乘值

%======== Here we are gonna chose the type of the incident light =========%在这里我们要选择入射光的类型
%=========================================================================%
if IB=="plane wave"            % Incident beam is plane wave平面波入射
    expikr = exp(1i*kr); %传播相位项
    Ex = E0(1)*expikr.*INDEX_INSIDE;
    Ey = E0(2)*expikr.*INDEX_INSIDE;
    Ez = E0(3)*expikr.*INDEX_INSIDE;

    E_x = reshape(Ex,[Nx,Ny,Nz]);
    E_y = reshape(Ey,[Nx,Ny,Nz]);
    E_z = reshape(Ez,[Nx,Ny,Nz]);
    E_vector = [Ex;Ey;Ez];   % (3N,1)

elseif IB=="gaussian"        % Incident beam is Gaussian 
    % 坐标选择：根据传播方向确定 r, z
    if K0(1)==1
        r = sqrt(r_block(:,2).^2 + r_block(:,3).^2);
        z = r_block(:,1);
    elseif K0(2)==1
        r = sqrt(r_block(:,1).^2 + r_block(:,3).^2);
        z = r_block(:,2);
    else
        r = sqrt(r_block(:,1).^2 + r_block(:,2).^2);
        z = r_block(:,3);
    end
    z=-z;
    W0 = Waist_r * Lambda;                    % 束腰半径
    K = 2*pi*nb / Lambda;                     % 波数
    zR = pi * (W0.^2) * nb / Lambda;          % 瑞利长度
    Qz = atan((z - z0)./zR);                  % Gouy 相位
    Wz = W0 .* sqrt(1 + ((z - z0)./zR).^2);   % 光斑大小参数
    Rz_inverse = z ./ ((z - z0).^2 + zR.^2);  % 曲率半径的倒数
 
    % 基模高斯场分布
    E = (W0./Wz).*exp(-(r./Wz).^2 - 1i*(K.*(z-z0) + K.*(r.^2).*Rz_inverse/2 - Qz));
    % E = (W0./Wz).*exp(-(r./Wz).^2 - 1i*(K.*(z-z0) + K.*(r.^2).*Rz_inverse/2 - Qz));
    Ex = E0(1)*E.*INDEX_INSIDE;
    Ey = E0(2)*E.*INDEX_INSIDE;
    Ez = E0(3)*E.*INDEX_INSIDE;

    E_x = reshape(Ex,[Nx,Ny,Nz]);
    E_y = reshape(Ey,[Nx,Ny,Nz]);
    E_z = reshape(Ez,[Nx,Ny,Nz]);
    E_vector = [Ex;Ey;Ez];

elseif IB=="laguerre-gaussian"   % Incident beam is Laguerre-Gaussian beam

    % -------- 坐标系选择（保持原逻辑） --------
    if K0(1)==1
        r = sqrt(r_block(:,2).^2 + r_block(:,3).^2);
        phi = atan2(r_block(:,3), r_block(:,2));
        z = r_block(:,1);
    elseif K0(2)==1
        r = sqrt(r_block(:,1).^2 + r_block(:,3).^2);
        phi = atan2(r_block(:,3), r_block(:,1));
        z = r_block(:,2);
    else
        r = sqrt(r_block(:,1).^2 + r_block(:,2).^2);
        phi = atan2(r_block(:,2), r_block(:,1));
        z = r_block(:,3);
    end
    z=-z;
    W0 = Waist_r * Lambda;              % 束腰半径
    k  = 2*pi*nb / Lambda;              % 波数
    zR = pi * W0^2 * nb / Lambda;       % Rayleigh length

    dz = z - z0;                        % 相对焦点坐标

    Wz = W0 .* sqrt(1 + (dz./zR).^2);   % 光斑大小
    Rz_inverse = dz ./ (dz.^2 + zR.^2); % 1 / R(z)
    Qz = atan(dz./zR);                  % Gouy 相位

    % -------- 归一化径向坐标 --------
    rho   = sqrt(2) .* r ./ Wz;
    alpha = abs(l);
    xarg  = rho.^2;

    % -------- 计算关联拉盖尔多项式 L_p^{|l|}(x) --------
    Lpl = zeros(size(xarg));
    for m = 0:p
        coef = (-1)^m * nchoosek(p + alpha, p - m) / factorial(m);
        Lpl = Lpl + coef .* (xarg.^m);
    end

    % -------- 轴附近数值处理 --------
    eps_r = max(W0*1e-8, 1e-12);
    phi(~isfinite(phi)) = 0;
    phi(r < eps_r) = 0;
    rho(r < eps_r) = 0;

    % -------- LG 光束复振幅（Gouy 相位符号已修正） --------
    norm_factor = sqrt(2*factorial(p)/(pi*factorial(p+abs(l))));
    amp_prefac = norm_factor * (1 ./ Wz) .* (rho.^alpha) .* Lpl .* exp(-rho.^2/2);
    phase = exp(-1i*(k.*dz + k.*(r.^2).*Rz_inverse/2) + 1i*(2*p+alpha+1).*Qz + 1i*l.*phi);
    E = amp_prefac .* phase;

    % -------- 不做任何额外归一化 --------
    Ex = E0(1) * E .* INDEX_INSIDE;
    Ey = E0(2) * E .* INDEX_INSIDE;
    Ez = E0(3) * E .* INDEX_INSIDE;

    E_x = reshape(Ex,[Nx,Ny,Nz]);
    E_y = reshape(Ey,[Nx,Ny,Nz]);
    E_z = reshape(Ez,[Nx,Ny,Nz]);
    E_vector = [Ex; Ey; Ez];

    % fprintf('LG(p=%d,l=%d) Gaussian-aligned: max|E|=%.3e, min|E|=%.3e\n', ...
    %     p, l, max(abs(E(:))), min(abs(E(:))));
end


end