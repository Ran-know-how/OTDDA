function [Ex_incident,Ey_incident,Ez_incident,Hx_incident,Hy_incident,Hz_incident]=Incident_Field_plane(r_plane,E0,H0,K0,z0,kvec,Waist_r,IB,Lambda,nb,p,l)
if IB=="plane wave"    % Incident beam is plane wave
    kr=kvec(1)*r_plane(:,1)+kvec(2)*r_plane(:,2)+kvec(3)*r_plane(:,3);
    expikr=exp(1i*kr);
    
    Ex_incident=(E0(1)*expikr);
    Ey_incident=(E0(2)*expikr);
    Ez_incident=(E0(3)*expikr);
    Hx_incident=(H0(1)*expikr);
    Hy_incident=(H0(2)*expikr);
    Hz_incident=(H0(3)*expikr);
elseif IB=="gaussian"   % Incident beam is Gaussian
    
    if K0(1)==1
        r=(r_plane(:,2).^2+r_plane(:,3).^2).^(0.5);
        z=r_plane(:,1);
    elseif K0(2)==1
        r=(r_plane(:,1).^2+r_plane(:,3).^2).^(0.5);
        z=r_plane(:,2);
    else
        r=(r_plane(:,1).^2+r_plane(:,2).^2).^(0.5);
        z=r_plane(:,3);
    end
    z=-z;
    W0=Waist_r*Lambda;                    %Waist radius of the beam
    K=2*pi*nb/Lambda;                     %wave number
    zR=pi*(W0.^2).*nb./Lambda;
    Qz=atan((z-z0)./zR);                    % Is the Gouy phase at z
    Wz=W0.*((1+((z-z0)./zR).^2).^0.5);      % the spot size parameter
    Rz_inverse=z./((z-z0).^2+zR.^2);        % inverse of radius of curvature
    
    E=(W0./Wz).*exp(-(r./Wz).^2-1i*(K.*(z-z0)+K.*(r.^2).*Rz_inverse/2-Qz));
    H=(W0./Wz).*exp(-(r./Wz).^2-1i*(K.*(z-z0)+K.*(r.^2).*Rz_inverse/2-Qz));
    Ex_incident=E0(1)*E;
    Ey_incident=E0(2)*E;
    Ez_incident=E0(3)*E;
    Hx_incident=H0(1)*H;
    Hy_incident=H0(2)*H;
    Hz_incident=H0(3)*H;
elseif IB=="laguerre-gaussian"   % Incident beam is Laguerre-Gaussian beam

    if K0(1)==1
        xcoord = r_plane(:,2);
        ycoord = r_plane(:,3);
        r   = sqrt(xcoord.^2 + ycoord.^2);
        phi = atan2(ycoord, xcoord);
        z   = r_plane(:,1);
    elseif K0(2)==1
        xcoord = r_plane(:,1);
        ycoord = r_plane(:,3);
        r   = sqrt(xcoord.^2 + ycoord.^2);
        phi = atan2(ycoord, xcoord);
        z   = r_plane(:,2);
    else
        xcoord = r_plane(:,1);
        ycoord = r_plane(:,2);
        r   = sqrt(xcoord.^2 + ycoord.^2);
        phi = atan2(ycoord, xcoord);
        z   = r_plane(:,3);
    end
    z=-z;
    % -------- 光束参数 --------
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
    % amp_prefac = (W0 ./ Wz) .* (rho.^alpha) .* Lpl .* exp(-rho.^2/2);
    phase = exp(-1i*(k.*dz + k.*(r.^2).*Rz_inverse/2) + 1i*(2*p+alpha+1).*Qz + 1i*l.*phi);
    % phase = exp(-1i*( ...
    %             k.*dz ...                         % 传播相位
    %           + k.*(r.^2).*Rz_inverse/2 ...       % 曲率相位
    %           + (2*p + alpha + 1).*Qz ...         % Gouy 相位（注意正号）
    %           + l.*phi ));                        % OAM 相位

    E = amp_prefac .* phase;
    % -------- 电场（入射） --------
    Ex_incident = E0(1) * E;
    Ey_incident = E0(2) * E;
    Ez_incident = E0(3) * E;

    % -------- 磁场（入射）：H = (k̂ × E) / η --------
    mu0  = 4*pi*1e-7;
    c0   = 299792458;
    eta_b = mu0 * c0 / nb;

    Hx_incident = ( K0(2).*Ez_incident - K0(3).*Ey_incident ) / eta_b;
    Hy_incident = ( K0(3).*Ex_incident - K0(1).*Ez_incident ) / eta_b;
    Hz_incident = ( K0(1).*Ey_incident - K0(2).*Ex_incident ) / eta_b;

end


end