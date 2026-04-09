function T = heat_solver_general(INDEX_IN, q_local, dx, k_p, h, T_inf)
% -------------------------------------------------------------------------
% 高性能版本：通用 3D 颗粒热求解器（无球体假设）
%   • FVM + Robin boundary
%   • 稀疏矩阵一次性写入
%   • 适合任意 shape（ellipsoid, rod, cube, NP cluster…）
% -------------------------------------------------------------------------

[Nx,Ny,Nz] = size(INDEX_IN);
idx3 = find(INDEX_IN);      % inside voxel list
N = length(idx3);

% 3D→1D mapping
map = zeros(Nx,Ny,Nz);
map(idx3) = 1:N;

% Laplacian conductance
alpha = k_p / dx^2;

% Directions
dirs = int32([1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1]);

% Preallocate triplets for sparse matrix
I = zeros(7*N,1);
J = zeros(7*N,1);
S = zeros(7*N,1);
b = zeros(N,1);

ptr = 1;

for n = 1:N
    % voxel coordinates
    [i,j,k] = ind2sub([Nx Ny Nz], idx3(n));

    diag_val = 0;
    rhs = q_local(i,j,k);

    for d = 1:6
        ni = i + dirs(d,1);
        nj = j + dirs(d,2);
        nk = k + dirs(d,3);

        % out of grid → convection boundary
        if ni<1 || ni>Nx || nj<1 || nj>Ny || nk<1 || nk>Nz
            diag_val = diag_val + (h/dx);
            rhs      = rhs + (h/dx)*T_inf;
            continue;
        end
        
        if INDEX_IN(ni,nj,nk)
            % conduction neighbor
            nid = map(ni,nj,nk);
            I(ptr) = n;
            J(ptr) = nid;
            S(ptr) = -alpha;
            ptr = ptr+1;

            diag_val = diag_val + alpha;
        else
            % convection boundary
            diag_val = diag_val + (h/dx);
            rhs      = rhs + (h/dx)*T_inf;
        end
    end

    % diagonal
    I(ptr) = n;
    J(ptr) = n;
    S(ptr) = diag_val;
    ptr = ptr+1;

    b(n) = rhs;
end

% trim unused tail
I = I(1:ptr-1);
J = J(1:ptr-1);
S = S(1:ptr-1);

A = sparse(I,J,S,N,N);

% Solve
Tvec = A\b;

% reconstruct 3D temperature
T = T_inf*ones(Nx,Ny,Nz);
T(idx3) = Tvec;
end
