function [Hx_out, Hy_out, Hz_out] = H_out(r_block, r_plane, N_plane, px, py, pz, kvec, Hx_incident, Hy_incident, Hz_incident)
    % 分块处理以减少内存占用，修复维度不匹配问题
    % 输入参数:
    %   r_block: [N^3, 3] 矩阵
    %   r_plane: [N^2, 3] 矩阵 (N_plane = N^2)
    %   其他参数保持不变

    % 计算k向量的模
    k_norm = norm(kvec);
    c = 3e8;

    % 初始化输出磁场（累加器）
    Hx = zeros(N_plane, 1);
    Hy = zeros(N_plane, 1);
    Hz = zeros(N_plane, 1);

    % 设置分块大小（根据内存情况调整）
    block_size = 1000;  % 可根据内存大小调整
    total_blocks = ceil(size(r_block, 1) / block_size);

    % 分块处理每个r_block子集
    for b = 1:total_blocks
        % 计算当前块的索引
        start_idx = (b-1)*block_size + 1;
        end_idx = min(b*block_size, size(r_block, 1));
        block_indices = start_idx:end_idx;
        current_block_size = end_idx - start_idx + 1;  % 当前块的实际大小

        % 提取当前块的数据
        r_block_sub = r_block(block_indices, :);
        px_sub = px(block_indices);
        py_sub = py(block_indices);
        pz_sub = pz(block_indices);

        % 调整维度以匹配广播要求
        % 将 [current_block_size, 1] 转换为 [1, current_block_size]
        px_sub = reshape(px_sub, 1, current_block_size);
        py_sub = reshape(py_sub, 1, current_block_size);
        pz_sub = reshape(pz_sub, 1, current_block_size);

        % 维度扩展
        r_plane_expanded = permute(r_plane, [1, 3, 2]);  % [N_plane, 1, 3]
        r_block_expanded = permute(r_block_sub, [3, 1, 2]);  % [1, current_block_size, 3]

        % 计算向量差 r = r_plane - r_block_sub [N_plane, current_block_size, 3]
        r = r_plane_expanded - r_block_expanded;

        % 计算r的长度 [N_plane, current_block_size]
        r_length = sqrt(sum(r.^2, 3));

        % 计算单位向量 [N_plane, current_block_size, 3]
        r_hat = r ./ r_length(:, :, ones(1, 3));  % 扩展维度匹配

        % 提取单位向量分量 [N_plane, current_block_size]
        rx = r_hat(:, :, 1);
        ry = r_hat(:, :, 2);
        rz = r_hat(:, :, 3);

        % 计算指数项和系数（当前块）
        Exp_ikvec_rjk = exp(1i * k_norm * r_length) ./ r_length;
        H_xishu = (c / (4 * pi)) * (1 - 1 ./ (1i * k_norm * r_length)) .* Exp_ikvec_rjk * k_norm^2;

        % 计算当前块对磁场的贡献（修复维度不匹配问题）
        Hx_contrib = sum((ry .* pz_sub - rz .* py_sub) .* H_xishu, 2);
        Hy_contrib = sum((rz .* px_sub - rx .* pz_sub) .* H_xishu, 2);
        Hz_contrib = sum((rx .* py_sub - ry .* px_sub) .* H_xishu, 2);

        % 累加贡献 
        Hx = Hx + Hx_contrib;
        Hy = Hy + Hy_contrib;
        Hz = Hz + Hz_contrib;

        % 清理内存
        clear r r_length r_hat Exp_ikvec_rjk H_xishu
        clear Hx_contrib Hy_contrib Hz_contrib
    end

    % 加上入射场
    % Hx_out = Hx + Hx_incident;
    % Hy_out = Hy + Hy_incident;
    % Hz_out = Hz + Hz_incident;
    Hx_out = Hx;
    Hy_out = Hy;
    Hz_out = Hz;
end