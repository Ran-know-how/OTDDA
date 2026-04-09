function [surf_idx, nx, ny, nz] = surface_normals_3D(INDEX_IN, dx)
% 计算三维表面法向量
% 输入:
%   INDEX_IN: 3D逻辑数组（内部为true，外部为false）
%   dx: 网格步长
% 输出:
%   surf_idx: 表面点线性索引
%   nx, ny, nz: 表面法向量分量（指向外部）

% ---- 确定GPU使用 ----
if isa(INDEX_IN,'gpuArray')
    toGPU = true;
    BW = gpuArray(double(INDEX_IN > 0));  % 确保BW在GPU上
else
    toGPU = false;
    BW = double(INDEX_IN > 0);
end
BW = gather(INDEX_IN);
% ---- 提取表面体素（边界） ----
if toGPU
    % 对于GPU，使用3D形态学操作提取表面
    % 创建3D结构元素
    se = strel3d(1);  % 3x3x3立方体结构元素
    BW_eroded = imerode(gather(BW), se);
    if toGPU
        BW_eroded = gpuArray(BW_eroded);
    end
    surf = BW & ~BW_eroded;
else
    % 对于CPU，使用bwperim
    surf = bwperim(BW, 26);  % 3D 26-connected boundary
end

% ---- 方法1：使用符号距离函数（SDF）的梯度（推荐） ----
% 计算符号距离函数：内部为负，外部为正
SDF = bwdist(~BW) - bwdist(BW);

% 如果使用GPU且需要转回CPU计算梯度
if toGPU
    SDF_cpu = gather(SDF);
else
    SDF_cpu = SDF;
end

% 计算SDF的梯度（MATLAB的gradient函数对3D数组的输出顺序是[gy, gx, gz]）
[gy, gx, gz] = gradient(SDF_cpu, dx, dx, dx);  % 注意：第一个输出是y方向梯度

% 归一化得到法向量（SDF梯度指向外部）
Gnorm = sqrt(gx.^2 + gy.^2 + gz.^2) + 1e-12;
nx = gx ./ Gnorm;  % x方向法向量分量
ny = gy ./ Gnorm;  % y方向法向量分量
nz = gz ./ Gnorm;  % z方向法向量分量

% ---- 方法2：使用二值图像梯度（备选，需要检查方向） ----
% 如果使用方法1有问题，可以尝试方法2，但需要验证方向
% 注意：BW的梯度方向可能指向内部，可能需要取反

% % 3D Sobel kernels
% Sx = cat(3, ...
%    [ -1 0 1; -3 0 3; -1 0 1 ], ...
%    [ -3 0 3; -6 0 6; -3 0 3 ], ...
%    [ -1 0 1; -3 0 3; -1 0 1 ]);
% 
% Sy = cat(3, ...
%    [ -1 -3 -1; 0 0 0; 1 3 1 ], ...
%    [ -3 -6 -3; 0 0 0; 3 6 3 ], ...
%    [ -1 -3 -1; 0 0 0; 1 3 1 ]);
% 
% Sz = cat(3, ...
%    [ -1 -3 -1; -3 -6 -3; -1 -3 -1 ], ...
%    [  0  0  0;  0  0  0;  0  0  0 ], ...
%    [  1  3  1;  3  6  3;  1  3  1 ]);
% 
% % 计算BW的梯度
% if toGPU
%     Sx = gpuArray(Sx); Sy = gpuArray(Sy); Sz = gpuArray(Sz);
% end
% 
% GX = convn(BW, Sx, 'same');
% GY = convn(BW, Sy, 'same');
% GZ = convn(BW, Sz, 'same');
% 
% % 归一化
% Gnorm = sqrt(GX.^2 + GY.^2 + GZ.^2) + 1e-12;
% nx_temp = -GX ./ Gnorm;  % 取负号，因为BW梯度指向内部
% ny_temp = -GY ./ Gnorm;
% nz_temp = -GZ ./ Gnorm;
% 
% % 验证方向：在凸区域，法向量应指向外部
% % 这里选择使用方法1的结果，但保留方法2作为参考

% ---- 提取表面点的法向量 ----
surf_idx = find(surf);

% 确保nx, ny, nz与surf_idx在同一设备上
if toGPU
    nx = gpuArray(nx);
    ny = gpuArray(ny);
    nz = gpuArray(nz);
end

nx = nx(surf_idx);
ny = ny(surf_idx);
nz = nz(surf_idx);

end

% ---- 辅助函数：创建3D结构元素 ----
function se = strel3d(r)
% 创建半径为r的3D球形结构元素（近似）
[x, y, z] = meshgrid(-r:r, -r:r, -r:r);
se = strel(sqrt(x.^2 + y.^2 + z.^2) <= r);
end