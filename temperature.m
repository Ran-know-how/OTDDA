function [T_sol,E_whole,inside_idx,XX,YY,ZZ] = temperature(r_eff,grid,Lambda,Np_shape,Structure,IB,E01,H01,K01,CT,GPU,np,nb)
eps0 = 8.85e-12;
mu0 = 4*pi*1e-7;
%处理相对介电常数
eps_nps=(np).^2;%粒子的
epsb=nb^2;%环境
if GPU==1
    eps_NP_eb=gpuArray(eps_nps./epsb);  % Ratio of metal-to-medium dielectric function
elseif GPU==0
    eps_NP_eb=eps_nps./epsb;  % Ratio of metal-to-medium dielectric function
end
%直径
d_eff=2*r_eff;
%最小分割单元的尺寸
d=r_eff/grid;
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
c0 = 3e8;
k=2*pi./Lambda*nb;
kvec=k*K0;
omega=2*pi*c0/Lambda;
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
%入射光场处理
if IB=="plane wave" 
    z0=0;              % Focus point of the Gaussian beam
    Waist_r=100;       % ratio of waist raduis of beam to wavelength
    p = 0;
    l= 0;
elseif IB=="gaussian"
    z0=0;              % Focus point of the Gaussian beam
    NA = 0.3;
    n_object = 1;
    Waist_r = n_object/(pi*NA);
    p = 0;
    l= 0;
elseif IB=="laguerre-gaussian"
    z0=0;              % Focus point of the Gaussian beam
    NA = 0.1;
    n_object = 1;
    Waist_r = n_object/(pi*NA);     
    p = 3;
    l= 0;
end 
%找到包络粒子的矩形
[Lx,Ly,Lz] = Nps_parameters(r_eff,Np_shape);
%xyz方向最大的长度，各个方向上偶极子的数量和总的偶极子的数量
[Max_x,Max_y,Max_z,N,Nx,Ny,Nz,r_block,X,Y,Z,d_inter]=Coordinates(GPU,d,Lx,Ly,Lz,...
    d_eff,Structure,arrangement); 
%给出偶极子的位置索引
[INDEX_INSIDE]=INDEX_INSIDE_NP(GPU,X,Y,Z,N,Np_shape,Lx,Ly,Lz,Structure,d_inter,arrangement);
INDEX_IN=reshape(INDEX_INSIDE,[Nx,Ny,Nz]);
%计算相互作用矩阵函数的组成部分
[rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,rjkrjk31_I,rjkrjk32_I,...
    rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,RJK]=RijRij(r_block);
%计算极化率（的倒数/逆）
[Inverse_Alpha]=Polarizability(GPU,kvec,eps_NP_eb,INDEX_IN,d,E0);
%给出入射电场
[E_x,E_y,E_z,E_vector]=Incident_Field(Lambda,IB,nb,r_block,kvec,K0,INDEX_INSIDE,Nx,Ny,Nz,E0,z0,Waist_r,p,l);
%相互作用矩阵表达式中的系数
Exp_ikvec_rjk=exp(1i*norm(kvec)*RJK)./RJK;
ikvec_rjk=(1i*norm(kvec)*RJK-1)./(RJK.^2);
%计算相互作用矩阵
[Axx,Axy,Axz,Ayy,Ayz,Azz]=Interaction_Matrix(kvec,Exp_ikvec_rjk,...
    ikvec_rjk, rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,...
    rjkrjk31_I,rjkrjk32_I,rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,Nx,Ny,Nz);
%傅里叶变换
[FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ]=FFT_Interaction(GPU,Axx...
    ,Axy,Axz,Ayy,Ayz,Azz,Nx,Ny,Nz);
%共轭梯度下降，得到偶极子强度
[px,py,pz]=Biconjugate_Gradient(E_x,E_y,E_z,Nx,Ny,Nz,N,Inverse_Alpha,...
    INDEX_IN,E_vector,FFT_AXX,FFT_AXY,FFT_AXZ,FFT_AYY,FFT_AYZ,FFT_AZZ,CT);
px=px.*INDEX_IN;
py=py.*INDEX_IN;
pz=pz.*INDEX_IN;
Ex_whole = Inverse_Alpha.*px;
Ey_whole = Inverse_Alpha.*py;
Ez_whole = Inverse_Alpha.*pz;

E_whole = sqrt((Ex_whole).*conj(Ex_whole)+(Ey_whole).*conj(Ey_whole)+Ez_whole.*conj(Ez_whole));

epsilon_2 = imag(eps_nps);  % 这里的 eps_nps 已经定义为复介电函数

q_local = 0.5 * omega .* epsilon_2 .* (abs(E_whole).^2) * eps0;

a = r_eff;   % simply rename for consistency

% ---------------- thermal/gas inputs (与第二段完全一致) ----------------
k_p = 0.2;
% k_p = 20;
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
% if ~isequal(size(q_local), size(XX))
%     q_local = reshape(q_local, size(XX));
% end

%% ------------------ 第二段中的球内 mask ----------------------------
mask = (r <= a);
inside_idx = find(mask);
Mvar = numel(inside_idx);
idmap = zeros(N*N*N,1,'int32');
idmap(inside_idx) = 1:Mvar;

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
    g = inside_idx(n);
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
T_sol(inside_idx) = T_vec;
end

