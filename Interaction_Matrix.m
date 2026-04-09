%=========================================================================%
%===== Calculating six tensor blocks: Axx, Axy, Axz, Ayy, Ayz and Azz ====%
%计算相互作用的六个分量
%=========================================================================%

function [Axx,Axy,Axz,Ayy,Ayz,Azz]=Interaction_Matrix(kvec,Exp_ikvec_rjk,...
    ikvec_rjk, rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,...
    rjkrjk31_I,rjkrjk32_I,rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,Nx,Ny,Nz)
epsilon0 = 8.85e-12;
xishu = 1/(4*pi*epsilon0);
A1=xishu*(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk1_I+ ikvec_rjk.*rjkrjk31_I));
% A1=(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk1_I+ ikvec_rjk.*rjkrjk31_I));
Axx=reshape(A1,[Nx,Ny,Nz]);
Axx(1,1,1)=0;
clear A1

A2=xishu*(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk2_I+ ikvec_rjk.*rjkrjk32_I));
% A2=(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk2_I+ ikvec_rjk.*rjkrjk32_I));
Axy=reshape(A2,[Nx,Ny,Nz]);
Axy(1,1,1)=0;
clear A2

A3=xishu*(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk3_I+ ikvec_rjk.*rjkrjk33_I));
% A3=(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk3_I+ ikvec_rjk.*rjkrjk33_I));
Axz=reshape(A3,[Nx,Ny,Nz]);
Axz(1,1,1)=0;
clear A3

A4=xishu*(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk4_I+ ikvec_rjk.*rjkrjk34_I));
% A4=(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk4_I+ ikvec_rjk.*rjkrjk34_I));
Ayy=reshape(A4,[Nx,Ny,Nz]);
Ayy(1,1,1)=0;
clear A4

A5=xishu*(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk5_I+ ikvec_rjk.*rjkrjk35_I));
% A5=(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk5_I+ ikvec_rjk.*rjkrjk35_I));
Ayz=reshape(A5,[Nx,Ny,Nz]);
Ayz(1,1,1)=0;
clear A5

A6=xishu*(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk6_I+ ikvec_rjk.*rjkrjk36_I));
% A6=(Exp_ikvec_rjk.*((norm(kvec)^2)*rjkrjk6_I+ ikvec_rjk.*rjkrjk36_I));
Azz=reshape(A6,[Nx,Ny,Nz]);
Azz(1,1,1)=0;
clear A6  Exp_ikvec_rjk

end