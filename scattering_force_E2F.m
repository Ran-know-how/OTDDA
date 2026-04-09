function F_vec = scattering_force_E2F(Ex_whole,Ey_whole,Ez_whole,Lambda,INDEX_IN,d)
    c0 = 3e8;
    omega=2*pi*c0./Lambda;
    eps0 = 8.85e-12;
    mu0 = 4*pi*1e-7;
    phase_Ex = unwrap(angle(Ex_whole), [], 1);  % 沿x方向解缠相位
    phase_Ey = unwrap(angle(Ey_whole), [], 2);  % 沿y方向
    phase_Ez = unwrap(angle(Ez_whole), [], 3);  % 沿z方向

    % 计算相位梯度（中心差分）
    phase_Ex(INDEX_IN == 0) = 0;
    phase_Ey(INDEX_IN == 0) = 0;
    phase_Ez(INDEX_IN == 0) = 0;
    [kx, ky, kz] = gradient_phase(phase_Ex, phase_Ey, phase_Ez, d, d, d);

    % 归一化波矢
    k_norm = sqrt(kx.^2 + ky.^2 + kz.^2);
    kx_hat = kx ./ (k_norm + eps); % 避免除以零
    ky_hat = ky ./ (k_norm + eps);
    kz_hat = kz ./ (k_norm + eps);
    k = omega / c0;  % 实际波数（c0为真空中光速）
    kx_hat = kx_hat * k;  
    ky_hat = ky_hat * k;
    kz_hat = kz_hat * k;

    E_field = cat(4, Ex_whole, Ey_whole, Ez_whole);
    k_field = cat(4, kx_hat, ky_hat, kz_hat);
    % 计算磁场（考虑局部波矢）
    Hx_vec = (k_field(:,:,:,2) .* E_field(:,:,:,3) - k_field(:,:,:,3) .* E_field(:,:,:,2)) / (omega * mu0).* INDEX_IN;
    Hy_vec = (k_field(:,:,:,3) .* E_field(:,:,:,1) - k_field(:,:,:,1) .* E_field(:,:,:,3)) / (omega * mu0).* INDEX_IN;
    Hz_vec = (k_field(:,:,:,1) .* E_field(:,:,:,2) - k_field(:,:,:,2) .* E_field(:,:,:,1)) / (omega * mu0).* INDEX_IN;
    H_field = cat(4, Hx_vec, Hy_vec, Hz_vec);
    F_vec = compute_MST_force(E_field, H_field, INDEX_IN, d);
end
