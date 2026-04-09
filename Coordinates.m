% Finding coordinate of the dipoles inside the rectangular block %求矩形块内偶极子的坐标
function [Max_x,Max_y,Max_z,N,Nx,Ny,Nz,r_block,X,Y,Z,d_inter]=Coordinates(GPU,d,Lx,Ly,Lz,...
                                                         d_eff,Structure,arrangement)
%输入量：d为单位体素大小，GPU为是否使用GPU,LxLyLz是包含NPs的矩形块尺寸，d_eff是有效尺寸，结构式单分子或者双分子，
dx=d;
dy=d;
dz=d;

if Structure=="monomeric" 
    d_inter=0;
    Max_x=Lx;
    Max_y=Ly;
    Max_z=Lz;
elseif Structure=="dimeric"           
    fprintf('\nA dimeric structure has been choosen.');
    d_ratio=input('\nEnter the ratio of (interparticle distance)/(effective diameter):');
    clc
    d_inter=d_ratio*d_eff;
                % For nanoparticles with head-tail orientation in z-direction
    if arrangement=="x_orient"  % side by side arrangement in x-direction
        Max_x=2*Lx+d_inter;
        Max_y=Ly;
        Max_z=Lz;
    elseif arrangement=="y_orient" % side by side arrangement in x-direction
        Max_x=Lx;
        Max_y=2*Ly+d_inter;
        Max_z=Lz;
    elseif arrangement=="z_orient"
        Max_x=Lx;    % head to tail arrangement in z-direction
        Max_y=Ly;
        Max_z=2*Lz+d_inter;
    end            
    
end


%======= Obtaining Coordinates of nanocubs or nanocells within NPS =======%
%=========================================================================%
if GPU==1
    ix=gpuArray(-round(Max_x/(2*dx)):round(Max_x/(2*dx)));
    iy=gpuArray(-round(Max_y/(2*dy)):round(Max_y/(2*dy)));
    iz=gpuArray(-round(Max_z/(2*dz)):round(Max_z/(2*dz)));
elseif GPU==0

    ix=-round(Max_x/(2*dx)):round(Max_x/(2*dx));
    iy=-round(Max_y/(2*dy)):round(Max_y/(2*dy));
    iz=-round(Max_z/(2*dz)):round(Max_z/(2*dz));
end

[y,x,z]=meshgrid(iy,ix,iz);
Nx=length(ix);                   % Number of the dipoles in the x-direction x方向的偶极子数量
Ny=length(iy);                   % Number of the dipoles in the y-direction y方向的偶极子数量
Nz=length(iz);                   % Number of the dipoles in the z-direction z方向的偶极子数量

N=Nx*Ny*Nz;%总的偶极子数量

X=reshape(x,[N,1])*dx; % X-coordinates of the dipoles 转化为一列便于计算
Y=reshape(y,[N,1])*dy; % Y-coordinates of the dipoles
Z=reshape(z,[N,1])*dz; % Z-coordinates of the dipoles

r_block=[X Y Z];         % Position of the each dipoles inside rectangular block 矩形块中间每一个偶极子的坐标
%=========================================================================%

end