function [Ex_out, Ey_out, Ez_out, Hx_out, Hy_out, Hz_out] = EM_out(r_block, r_plane, N_plane, px, py, pz, kvec, Ex_incident, Ey_incident, Ez_incident, Hx_incident, Hy_incident, Hz_incident)
    % 同时计算电场E与磁场H的离散偶极子近似（DDA）函数
    % 输入参数:
    %   r_block: [N^3, 3] 偶极子位置矩阵（N为粒子离散网格数）
    %   r_plane: [N^2, 3] 观测点位置矩阵（N_plane = N^2）
    %   px/py/pz: [N^3, 1] 偶极矩分量（与r_block一一对应）
    %   kvec: [3, 1] 波矢（如入射波方向）
    %   Ex_incident/Ey_incident/Ez_incident: [N^2, 1] 入射电场（需叠加入总电场）
    %   Hx_incident/Hy_incident/Hz_incident: [N^2, 1] 入射磁场（需叠加入总磁场）
    % 输出参数:
    %   Ex_out/Ey_out/Ez_out: [N^2, 1] 总电场分量（散射场+入射场）
    %   Hx_out/Hy_out/Hz_out: [N^2, 1] 总磁场分量（散射场+入射场）
    
    % ------------------------------
    % 1. 基础参数与初始化
    % ------------------------------
    k_norm = norm(kvec);                % 波数（k = |kvec|）
    eps0 = 8.854187817e-12;             % 真空介电常数（SI单位）
    c = 3e8;                            % 光速（SI单位）
    
    % 初始化电磁场累加器（散射场）
    Ex = zeros(N_plane, 1);             % 电场x分量（散射场）
    Ey = zeros(N_plane, 1);             % 电场y分量（散射场）
    Ez = zeros(N_plane, 1);             % 电场z分量（散射场）
    Hx = zeros(N_plane, 1);             % 磁场x分量（散射场）
    Hy = zeros(N_plane, 1);             % 磁场y分量（散射场）
    Hz = zeros(N_plane, 1);             % 磁场z分量（散射场）
    
    % ------------------------------
    % 2. 分块处理（减少内存占用）
    % ------------------------------
    block_size = 2000;                  % 每块偶极子数量（可根据内存调整）
    total_blocks = ceil(size(r_block, 1) / block_size);  % 总块数

    for b = 1:total_blocks
        % 2.1 提取当前块数据
        start_idx = (b-1)*block_size + 1;
        end_idx = min(b*block_size, size(r_block, 1));
        r_block_sub = r_block(start_idx:end_idx, :);  % 当前块偶极子位置
        
        % 调整偶极矩为行向量（广播计算用）
        px_sub = reshape(px(start_idx:end_idx), 1, []);
        py_sub = reshape(py(start_idx:end_idx), 1, []);
        pz_sub = reshape(pz(start_idx:end_idx), 1, []);
        M = size(px_sub, 2);             % 当前块偶极子数量
        
        % 2.2 维度扩展（广播计算向量差）
        r_plane_exp = permute(r_plane, [1, 3, 2]);  % [N_plane, 1, 3]（观测点扩展）
        r_block_exp = permute(r_block_sub, [3, 1, 2]);  % [1, M, 3]（偶极子扩展）
        r = r_plane_exp - r_block_exp;  % 向量差：r = 观测点 - 偶极子（[N_plane, M, 3]）
        
        
        % 2.3 计算距离（处理零值，避免除以零）
        r_length = sqrt(sum(r.^2, 3));  % 距离：R = |r|（[N_plane, M]）
        r_safe = max(r_length, eps);    % 替换零值为eps（避免除以零）
        % 2.4 计算单位向量（\hat{n} = r/R）
        r_hat = r ./ r_safe(:, :, ones(1, 3));  % 单位向量（[N_plane, M, 3]）
        rx = r_hat(:, :, 1);             % 单位向量x分量（[N_plane, M]）
        ry = r_hat(:, :, 2);             % 单位向量y分量（[N_plane, M]）
        rz = r_hat(:, :, 3);             % 单位向量z分量（[N_plane, M]）
        
        % 2.5 计算偶极矩与单位向量点积（\hat{n}·p）
        p_dot_r = rx .* px_sub + ry .* py_sub + rz .* pz_sub;  % [N_plane, M]
        
        % 2.6 计算相位因子（e^(ikR)）
        exp_ikr = exp(1i * k_norm * r_safe);  % 相位因子（[N_plane, M]）
        
        % ------------------------------
        % 3. 计算电场贡献（E）
        % ------------------------------
        % 3.1 远场项（k²项）：k² (p - \hat{n}(\hat{n}·p))（向量恒等式展开）
        term1_x = k_norm^2 * (px_sub - rx .* p_dot_r);  % 远场x分量（[N_plane, M]）
        term1_y = k_norm^2 * (py_sub - ry .* p_dot_r);  % 远场y分量（[N_plane, M]）
        term1_z = k_norm^2 * (pz_sub - rz .* p_dot_r);  % 远场z分量（[N_plane, M]）
        
        % 3.2 近场项（3\hat{n}(\hat{n}·p)项）：3(1 - jkR)/R² * \hat{n}(\hat{n}·p)
        term2_coeff = 3 * (1 - 1i * k_norm * r_safe) ./ (r_safe.^2);  % 近场系数（[N_plane, M]）
        term2_x = term2_coeff .* rx .* p_dot_r;  % 近场x分量（[N_plane, M]）
        term2_y = term2_coeff .* ry .* p_dot_r;  % 近场y分量（[N_plane, M]）
        term2_z = term2_coeff .* rz .* p_dot_r;  % 近场z分量（[N_plane, M]）
        
        % 3.3 合并远场与近场项
        total_E_x = term1_x + term2_x;  % 电场x分量总项（[N_plane, M]）
        total_E_y = term1_y + term2_y;  % 电场y分量总项（[N_plane, M]）
        total_E_z = term1_z + term2_z;  % 电场z分量总项（[N_plane, M]）
        
        % 3.4 电场公共系数（1/(4πε₀) * e^(ikR)/R）
        common_coeff_E = exp_ikr ./ (4 * pi * eps0 * r_safe);  % [N_plane, M]
        
        % 3.5 累加当前块电场贡献
        Ex = Ex + sum(total_E_x .* common_coeff_E, 2);
        Ey = Ey + sum(total_E_y .* common_coeff_E, 2);
        Ez = Ez + sum(total_E_z .* common_coeff_E, 2);
        
        % ------------------------------
        % 4. 计算磁场贡献（H）
        % ------------------------------
        % 4.1 磁场系数（对应原H_out函数的H_xishu）
        Exp_ikvec_rjk = exp_ikr ./ r_safe;  % 相位因子与距离比（[N_plane, M]）
        H_coeff = (c / (4 * pi)) * (1 - 1 ./ ( - 1i * k_norm * r_safe)) .* Exp_ikvec_rjk * k_norm^2;  % [N_plane, M]
        
        % 4.2 计算磁场交叉项（对应原H_out函数的(ry*pz - rz*py)等）
        cross_H_x = ry .* pz_sub - rz .* py_sub;  % 磁场x分量交叉项（[N_plane, M]）
        cross_H_y = rz .* px_sub - rx .* pz_sub;  % 磁场y分量交叉项（[N_plane, M]）
        cross_H_z = rx .* py_sub - ry .* px_sub;  % 磁场z分量交叉项（[N_plane, M]）
        
        % 4.3 累加当前块磁场贡献
        Hx = Hx + sum(cross_H_x .* H_coeff, 2);
        Hy = Hy + sum(cross_H_y .* H_coeff, 2);
        Hz = Hz + sum(cross_H_z .* H_coeff, 2);
        
        % ------------------------------
        % 5. 清理当前块中间变量（减少内存占用）
        % ------------------------------
        clear r_plane_exp r_block_exp r r_length r_safe r_hat rx ry rz p_dot_r exp_ikr
        clear term1_x term1_y term1_z term2_coeff term2_x term2_y term2_z total_E_x total_E_y total_E_z common_coeff_E
        clear Exp_ikvec_rjk H_coeff cross_H_x cross_H_y cross_H_z
    end
    
    % ------------------------------
    % 6. 叠加入射场（总场 = 散射场 + 入射场）
    % ------------------------------
    % Ex_out = Ex + Ex_incident;
    % Ey_out = Ey + Ey_incident;
    % Ez_out = Ez + Ez_incident;
    % Hx_out = Hx + Hx_incident;
    % Hy_out = Hy + Hy_incident;
    % Hz_out = Hz + Hz_incident;
    Ex_out = Ex;
    Ey_out = Ey;
    Ez_out = Ez;
    Hx_out = Hx;
    Hy_out = Hy;
    Hz_out = Hz;
end