function [Fx,Fy,Fz] = Photophoretic_force(T_sol,XX,YY,ZZ,r_eff,grid)
% % ------------------ Photophoretic force (approx) -----------------
% Assumptions: ambient pressure p_amb (Pa), ideal gas, diffuse re-emission
p_amb = 101325;            % ambient pressure, modify if needed
T_inf = 300;               % ambient temperature (K) from your earlier
kB = 1.380649e-23;
m_air = 4.81e-26;          % mean molecular mass air ~ 28.97 amu -> 28.97*1.6605e-27 ~ 4.81e-26 kg
lambda_mfp = 65e-9;   
d = r_eff/grid;
r = sqrt(XX.^2 + YY.^2 + ZZ.^2);
% surface sampling: we used surf_idx earlier (points with abs(r-a) <= surf_tol)
if ~exist('surf_idx','var') || isempty(surf_idx)
    surf_tol = d*1.2;
    surf_idx = find(abs(r - r_eff) <= surf_tol);
end
% positions and normals of surface sample points
xs = XX(surf_idx); ys = YY(surf_idx); zs = ZZ(surf_idx);
rvecs = [xs(:), ys(:), zs(:)];
rnorms = sqrt(sum(rvecs.^2,2));
nvecs = rvecs ./ rnorms;    % outward normal (approx)

% local surface temperature
Tsurf = T_sol(surf_idx);

% estimate area element for each surface sampling point
Nsurf = numel(surf_idx);
A_sphere = 4*pi*r_eff^2;
dA = A_sphere / Nsurf;    % uniform approx; improve by Voronoi on sphere if needed

% mean gas number density
n_g = p_amb / (kB * T_inf);

% average thermal velocity (not strictly needed for our simple model)
c_bar = sqrt(8*kB*T_inf/(pi*m_air));

% Knudsen number (you computed earlier)
Kn = lambda_mfp / r_eff;

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
end
