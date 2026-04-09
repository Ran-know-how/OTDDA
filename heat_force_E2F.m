function [F_vec,T_sol] = heat_force_E2F(r_eff,E_whole,eps_nps,Lambda,Nx,Ny,Nz)

c0 = 3e8;
omega=2*pi*c0/Lambda;
epsilon_2 = imag(eps_nps);  % 这里的 eps_nps 已经定义为复介电函数
eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;
q_local = 0.5 * omega .* epsilon_2 .* (abs(E_whole).^2) * eps0;

a = r_eff;   % simply rename for consistency

% ---------------- thermal/gas inputs (与第二段完全一致) ----------------
k_p = 0.2;            
k_g = 0.0257;         
lambda_mfp = 65e-9;   
beta = 0.8;           
T_inf = 300;          
Kn = lambda_mfp / a;
h_cont = k_g / a;
h_fm = k_g / (beta * lambda_mfp);
h_eff = (1 ./ (1./h_cont + 1./h_fm));

fprintf('a = %.1e m, Kn = %.3f\n', a, Kn);
fprintf('h_cont = %.3e, h_fm = %.3e, h_eff = %.3e\n', h_cont, h_fm, h_eff);

%% ------------------ 构建第二段中的构造网格 --------------------------
% 第二段用的是：grid = linspace(-a, a, N);
% 第一段使用 Nx = Ny = Nz
N = Nx;
L = a;
grid = linspace(-L, L, N);
dx = grid(2) - grid(1);
[XX,YY,ZZ] = ndgrid(grid,grid,grid);
r = sqrt(XX.^2 + YY.^2 + ZZ.^2);
% reshape q_local 以保证其 shape 一致
if ~isequal(size(q_local), size(XX))
    q_local = reshape(q_local, size(XX));
end

%% ------------------ 第二段中的球内 mask ----------------------------
mask = (r <= a);
unknown_idx = find(mask);
Mvar = numel(unknown_idx);
idmap = zeros(N*N*N,1,'int32');
idmap(unknown_idx) = 1:Mvar;

%% ------------------ 第二段的稀疏系数矩阵构建 ----------------------
coef_off = -k_p/(dx^2);
max_entries = Mvar * 7;

I = zeros(max_entries,1);
J = zeros(max_entries,1);
V = zeros(max_entries,1);
b = zeros(Mvar,1);
ptr = 0;

idx_lin = @(ii,jj,kk) sub2ind([N,N,N], ii,jj,kk);
eps_small = 1e-12;

for n = 1:Mvar
    g = unknown_idx(n);
    [ii,jj,kk] = ind2sub([N,N,N], g);
    xi = grid(ii); yi = grid(jj); zi = grid(kk);
    diag_val = 0;
    rhs = q_local(ii,jj,kk);

    neigh = [ii+1,jj,kk, 1,0,0;
             ii-1,jj,kk,-1,0,0;
             ii,jj+1,kk,0,1,0;
             ii,jj-1,kk,0,-1,0;
             ii,jj,kk+1,0,0,1;
             ii,jj,kk-1,0,0,-1];
    
    for m = 1:6
        ni = neigh(m,1); nj = neigh(m,2); nk = neigh(m,3);
        ux = neigh(m,4); uy = neigh(m,5); uz = neigh(m,6);

        outside = (ni < 1 || ni > N || nj < 1 || nj > N || nk < 1 || nk > N || ~mask(ni,nj,nk));

        if ~outside
            % internal neighbor
            ptr = ptr + 1;
            I(ptr)=n; J(ptr)=idmap(idx_lin(ni,nj,nk)); V(ptr)=coef_off;
            diag_val = diag_val + k_p/(dx^2);
        else
            % Robin boundary (第二段的严格写法)
            if ux~=0
                R = a^2 - yi^2 - zi^2;
                if R<0, t=dx; else, t = (ux>0)*(-xi+sqrt(R)) + (ux<0)*(xi+sqrt(R)); end
            elseif uy~=0
                R = a^2 - xi^2 - zi^2;
                if R<0, t=dx; else, t = (uy>0)*(-yi+sqrt(R)) + (uy<0)*(yi+sqrt(R)); end
            else
                R = a^2 - xi^2 - yi^2;
                if R<0, t=dx; else, t = (uz>0)*(-zi+sqrt(R)) + (uz<0)*(zi+sqrt(R)); end
            end
            
            d_i = max(eps_small, min(t, dx*1.5));

            face_coeff = (k_p * h_eff) / ( dx * (k_p + h_eff * d_i) );

            diag_val = diag_val + face_coeff;
            rhs = rhs + face_coeff * T_inf;
        end
    end

    ptr = ptr+1;
    I(ptr)=n; J(ptr)=n; V(ptr)=diag_val;
    b(n) = rhs;
end

I = I(1:ptr); J = J(1:ptr); V = V(1:ptr);
A = sparse(I,J,V,Mvar,Mvar);

%% ------------------- 求解温度场（第二段方法） ---------------------
if Mvar <= 6e4
    T_vec = A\b;
else
    tol = 1e-8; maxit = 2000;
    try
        L = ichol(A, struct('droptol',1e-3,'diagcomp',1e-3));
        [T_vec,flag] = pcg(A,b,tol,maxit,L,L');
    catch
        T_vec = pcg(A,b,tol,maxit);
    end
end

T_sol = ones(N,N,N)*T_inf;
T_sol(unknown_idx) = T_vec;
%画出了电场和温度的分布图
% figure('Name','Temperature and |E|^2 inside sphere','Position',[100 100 1200 500]);
% subplot(1,2,1);
% inside_idx = find(mask);
% scatter3(XX(inside_idx)*1e6, YY(inside_idx)*1e6, ZZ(inside_idx)*1e6, 20, T_sol(inside_idx), 'filled');
% axis equal; colorbar; title('T inside (K)'); xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)');
% view(45,25);
% 
% subplot(1,2,2);
% scatter3(XX(inside_idx)*1e6, YY(inside_idx)*1e6, ZZ(inside_idx)*1e6, 20, E_whole(inside_idx), 'filled');
% axis equal; colorbar; title('|E|^2 inside (V^2/m^2)'); xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)');
% view(45,25);

%% ------------------ 输出 ----------------------
fprintf('T_inside: min=%g, max=%g, mean=%g\n', ...
    min(T_sol(mask)), max(T_sol(mask)), mean(T_sol(mask)));
% % ------------------ Photophoretic force (approx) -----------------
% Assumptions: ambient pressure p_amb (Pa), ideal gas, diffuse re-emission
p_amb = 101325;            % ambient pressure, modify if needed
T_inf = 300;               % ambient temperature (K) from your earlier
kB = 1.380649e-23;
m_air = 4.81e-26;          % mean molecular mass air ~ 28.97 amu -> 28.97*1.6605e-27 ~ 4.81e-26 kg

% surface sampling: we used surf_idx earlier (points with abs(r-a) <= surf_tol)
if ~exist('surf_idx','var') || isempty(surf_idx)
    surf_tol = dx*1.2;
    surf_idx = find(abs(r - a) <= surf_tol);
end
% positions and normals of surface sample points
xs = XX(surf_idx); ys = YY(surf_idx); zs = ZZ(surf_idx);
rvecs = [xs(:), ys(:), zs(:)];
rnorms = sqrt(sum(rvecs.^2,2));
nvecs = rvecs ./ max(rnorms,eps);    % outward normal (approx)

% local surface temperature
Tsurf = T_sol(surf_idx);

% estimate area element for each surface sampling point
Nsurf = numel(surf_idx);
A_sphere = 4*pi*a^2;
dA = A_sphere / Nsurf;    % uniform approx; improve by Voronoi on sphere if needed

% mean gas number density
n_g = p_amb / (kB * T_inf);

% average thermal velocity (not strictly needed for our simple model)
c_bar = sqrt(8*kB*T_inf/(pi*m_air));

% Knudsen number (you computed earlier)
Kn = lambda_mfp / a;

% bridging / empirical coefficient: tuneable (default: smoothly reduces to 0 as Kn->infty)
% (This is a heuristic bridge — replace with Rohatschek expression for more accuracy.)
Cbridge = 1 ./ (1 + 10*Kn);   % << adjustable factor: 10 is empirical; smaller->closer to free-molecular

% local "effective pressure difference per area" model (engineering approx)
% deltaP_i = p_amb * (Tinf - Tsurf_i) / Tinf * Cbridge
deltaP = p_amb .* (T_inf - Tsurf) ./ T_inf .* Cbridge;

% compute force by summing vector contributions
Fx = sum( nvecs(:,1) .* deltaP * dA );
Fy = sum( nvecs(:,2) .* deltaP * dA );
Fz = sum( nvecs(:,3) .* deltaP * dA );
F_vec = [Fx, Fy, Fz];
F_vec = -1*F_vec;
end