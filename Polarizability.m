% %=========================================================================%
% %%Obtaining inverse of polarizibilty of each nanocube at different Lambda%
% %获得每个纳米立方体在不同波长下的极性斜率的倒数%
% %=========================================================================%
% 
function [Inverse_Alpha]=Polarizability(GPU,kvec,eps_NP_eb,INDEX_IN,d,E0)
%输入参量：GPU标志，波矢，金属与介质介电函数之比，偶极子索引，单位体素大小，电场振动方向
epsilon0 = 8.85e-12;
b1 = -1.891531;
b2 = 0.1648469;
b3 = -1.7700004;
dcube = d^3;
a_hat = kvec/norm(kvec);%norm,求范数,norm(m)是求欧几里得范数（二阶），及所有分量平方和。这里因此得到K的单位向量（方向）
e_hat = E0/norm(E0);%得到E的单位向量（反向）
S = 0;

for j = 1:3
    S = S + (a_hat(j)*e_hat(j))^2;%单位向量点乘
end
%S为一个数
if GPU==1
    S=gpuArray(S);
end

% a_CM = 3*dcube/(4*pi)*(eps_NP_eb - 1)./(eps_NP_eb + 2);%高斯制
a_CM =3*epsilon0*dcube*(eps_NP_eb - 1)./(eps_NP_eb + 2); % Clausius-Mossotti克劳修斯-莫索提关系
%见论文“MATLAB package for discrete dipole approximation by graphics
%processing unit: Fast Fourier Transform and Biconjugate Gradient”中的公式
%为一个复数
anr=a_CM./(1 + (a_CM/dcube).*(b1+eps_NP_eb*b2+eps_NP_eb*b3*S)*((norm(kvec)*d)^2)-(2/3)*1i*(norm(kvec)*dcube)^3);%晶格色散极化关系
%为一个复数

%S含义见晶格色散极化关系的公式

if GPU==1     % running in GPU
    aLDR=gpuArray(anr./(1-2/3*1i*(anr/dcube)*((norm(kvec)*d)^3)));
elseif GPU==0 % running in CPU
    aLDR=(anr./(1-2/3*1i*(anr/dcube)*((norm(kvec)*d)^3)));
end
% 修正的Lorentz-Lorenz极化率（含retardation）
% alpha = 3*epsilon0*dcube .* (eps_NP_eb - 1) ./ ...
%         (eps_NP_eb + 2 - (1i*norm(kvec)*d)^2/3 .* (eps_NP_eb - 1));
Inverse_Alpha=1/aLDR*INDEX_IN;%极化率的逆
% Inverse_Alpha=1/alpha*INDEX_IN;%极化率的逆
end
% function [Inverse_Alpha_test] = Polarizability(GPU,kvec,eps_NP_eb,INDEX_IN,d,E0)
% epsilon0 = 8.854e-12;
% dcube = d^3;
% a_CM = 3*epsilon0*dcube*(eps_NP_eb - 1)./(eps_NP_eb + 2);  % SI units
% 
% % 逐元素倒数
% Inverse_Alpha_test = INDEX_IN ./ a_CM;
% 
% end

