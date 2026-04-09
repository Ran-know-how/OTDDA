function [r_block,r_plane,Nx,Ny,Nz,N_plane,x_plane,y_plane,z_plane,SIZE,Xex,Yex,Zex]=...
    Two_D_plane_for_forwardscattering(GPU,d,Max_x,Max_y,Max_z,Lx,Ly,Lz,plane_2D)

dx=d;
dy=d;
dz=d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if contains(plane_2D,"xy")
    distance = Max_z;
    if GPU==1
        X1=gpuArray(-round(0.5*Max_x/dx+Lx/(2*dx)):round(0.5*Max_x/dx+Lx/(2*dx)));
        %[1，4Nx]的矩阵，坐标轴，取值范围是-2Nx到2Nx，得到了更大的网格
        X_extend=X1*dx;
        %乘上dx变成坐标
        Y1=gpuArray(-round(0.5*Max_y/dy+Ly/(2*dy)):round(0.5*Max_y/dy+Ly/(2*dy)));
        Y_extend=Y1*dy;
        Z1=gpuArray(-round(Max_z/(dz)):round(Max_z/(dz)));
        Z_extend=Z1*dz;
        %[1,4Nz]的矩阵，取值范围是-Nz到Nz，非作图平面维持原样

    elseif GPU==0
        X1=-round(0.5*Max_x/dx+Lx/(2*dx)):round(0.5*Max_x/dx+Lx/(2*dx));
        X_extend=X1*dx;
        Y1=-round(0.5*Max_y/dy+Ly/(2*dy)):round(0.5*Max_y/dy+Ly/(2*dy));
        Y_extend=Y1*dy;
        Z1=-round(Max_z/(dz)):round(Max_z/(dz));
        Z_extend=Z1*dz;
    end

    % [Yex,Xex]=meshgrid(Y_extend,X_extend);
    [Xex,Yex]=meshgrid(X_extend,Y_extend);
    SIZE=size(Xex);%（1，2）的double，里面是4Ny，4Nz
    if GPU==1
        Zex=zeros(SIZE(1,1),SIZE(1,2),'gpuArray');
    else
        Zex=zeros(SIZE(1,1),SIZE(1,2));
    end
    %在GPU中生成一个网格
    N_plane=SIZE(1,1)*SIZE(1,2);   % Total number of voxeles within NPs and in outside regionNP内和外部区域中的体素总数
    %目标平面上的点的总数
%     Nxex=length(X_extend);    % Number of voxels in the extended x-axis扩展x轴中的体素数量
%     Nyex=Ny;
%     Nzex=length(Z_extend);    % Number of voxels in the extended z-axis扩展z轴中的体素数量
    x_plane=reshape(Xex,[N_plane,1]);   % X coordinate of the voxels in the extended x-axis体素在扩展X轴上的X坐标
    y_plane=reshape(Yex,[N_plane,1]);   % y coordinate of the voxels in the extended y-axis体素在扩展Y轴上的X坐标
    z_plane=reshape(Zex,[N_plane,1]);   % z coordinate of the voxels in the extended z-axis体素在扩展Z轴上的X坐标
    % [y_rec,x_rec,z_rec]=meshgrid(Y_extend,X_extend,Z_extend);
    [x_rec,y_rec,z_rec]=meshgrid(X_extend,Y_extend,Z_extend);
   %======================================================================% 
elseif contains(plane_2D,"xz")
    distance = Max_y;
    if GPU==1
        X1=gpuArray(-round(0.5*Max_x/dx+Lx/(2*dx)):round(0.5*Max_x/dx+Lx/(2*dx)));
        X_extend=X1*dx;
        Y1=gpuArray(-round(Max_y/(dy)):round(Max_y/(dy)));
        Y_extend=Y1*dy;
        Z1=gpuArray(-round(0.5*Max_z/dz+Lz/(2*dz)):round(0.5*Max_z/dz+Lz/(2*dz)));
        Z_extend=Z1*dz;

    elseif GPU==0
        X1=-round(0.5*Max_x/dx+Lx/(2*dx)):round(0.5*Max_x/dx+Lx/(2*dx));
        X_extend=X1*dx;
        Y1=-round(Max_y/(dy)):round(Max_y/(dy));
        Y_extend=Y1*dy;
        Z1=-round(0.5*Max_z/dz+Lz/(2*dz)):round(0.5*Max_z/dz+Lz/(2*dz));
        Z_extend=Z1*dz;
    end
    % [Zex,Xex]=meshgrid(Z_extend,X_extend);
    [Xex,Zex]=meshgrid(X_extend,Z_extend);
    SIZE=size(Xex);
    if GPU==1
        Yex=zeros(SIZE(1,1),SIZE(1,2),'gpuArray');
    else
        Yex=zeros(SIZE(1,1),SIZE(1,2));
    end
    N_plane=SIZE(1,1)*SIZE(1,2);   % Total number of voxeles within NPs and in outside region

%     Nxex=length(X_extend);    % Number of voxels in the extended x-axis
%     Nyex=Ny;
%     Nzex=length(Z_extend);    % Number of voxels in the extended z-axis
    x_plane=reshape(Xex,[N_plane,1]);   % X coordinate of the voxels in the extended x-axis
    y_plane=reshape(Yex,[N_plane,1]);   % y coordinate of the voxels in the extended y-axis
    z_plane=reshape(Zex,[N_plane,1]);   % z coordinate of the voxels in the extended z-axis
    %放在一列方便运算
    % [y_rec,x_rec,z_rec]=meshgrid(Y_extend,X_extend,Z_extend);
    [x_rec,y_rec,z_rec]=meshgrid(X_extend,Y_extend,Z_extend);
    %拓展之后的体网格，格式为[4Ny,4Nx,2Nz]
    %=====================================================================%
elseif contains(plane_2D,"yz")
    distance = Max_z;
    if GPU==1
        X1=gpuArray(-round(Max_x/(dx)):round(Max_x/(dx))); 
        X_extend=X1*dx;
        Y1=gpuArray(-round(0.5*Max_y/dy+Ly/(2*dy)):round(0.5*Max_y/dy+Ly/(2*dy)));
        Y_extend=Y1*dy;
        Z1=gpuArray(-round(0.5*Max_z/dz+Lz/(2*dz)):round(0.5*Max_z/dz+Lz/(2*dz)));
        Z_extend=Z1*dz;

    elseif GPU==0
        X1=-round(Max_x/(dx)):round(Max_x/(dx)); 
        X_extend=X1*dx;
        Y1=-round(0.5*Max_y/dy+Ly/(2*dy)):round(0.5*Max_y/dy+Ly/(2*dy));
        Y_extend=Y1*dy;
        Z1=-round(0.5*Max_z/dz+Lz/(2*dz)):round(0.5*Max_z/dz+Lz/(2*dz));
        Z_extend=Z1*dz;
    end

    % [Zex,Yex]=meshgrid(Z_extend,Y_extend);
    [Yex,Zex]=meshgrid(Y_extend,Z_extend);
    SIZE=size(Yex);
    if GPU==1
        Xex=zeros(SIZE(1,1),SIZE(1,2),'gpuArray');
    else
        Xex=zeros(SIZE(1,1),SIZE(1,2));
    end
    N_plane=SIZE(1,1)*SIZE(1,2);   % Total number of voxeles within NPs and in outside regionNP内和外部区域中的体素总数
    %Nxex=length(X_extend);    % Number of voxels in the extended x-axis扩展x轴中的体素数量
    %Nyex=Ny;
    %Nzex=length(Z_extend);    % Number of voxels in the extended z-axis扩展z轴中的体素数量

    x_plane=reshape(Xex,[N_plane,1]);   % X coordinate of the voxels in the extended x-axis体素在扩展X轴上的X坐标
    y_plane=reshape(Yex,[N_plane,1]);   % y coordinate of the voxels in the extended y-axis体素在扩展Y轴上的X坐标
    z_plane=reshape(Zex,[N_plane,1]);   % z coordinate of the voxels in the extended z-axis体素在扩展Z轴上的X坐标
    % [y_rec,x_rec,z_rec]=meshgrid(Y_extend,X_extend,Z_extend);
    [x_rec,y_rec,z_rec]=meshgrid(X_extend,Y_extend,Z_extend);
end


r_plane=[x_plane y_plane z_plane];            % Position of the each nanocubes inside and outside of the NPs boundary
if plane_2D=="xy"
    r_plane(:,3) = 0;
elseif plane_2D=="xz"
    r_plane(:,2) = 0;
elseif plane_2D=="yz"
r_plane(:,1) = 0;
elseif plane_2D=="yz_x=2N"
r_plane(:,1) = distance;
x_plane = r_plane(:,1);
elseif plane_2D=="yz_x=-2N"
r_plane(:,1) = -distance;
x_plane =r_plane(:,1);
elseif plane_2D=="xy_z=2N"
r_plane(:,3) = distance;
z_plane = r_plane(:,3);
elseif plane_2D=="xy_z=-2N"
r_plane(:,3) = -distance;
z_plane= r_plane(:,3);
elseif plane_2D=="xz_y=2N"
r_plane(:,2) = distance;
y_plane= r_plane(:,2);
elseif plane_2D=="xz_y=-2N"
r_plane(:,2) = -distance;
y_plane= r_plane(:,2);
end
% [y_rec,x_rec,z_rec]=meshgrid(Y_extend,X_extend,Z_extend);
% subplot(131)
% plot(r_plane(:,1))
% subplot(132)
% plot(r_plane(:,2))
% subplot(133)
% plot(r_plane(:,3))

Nx=length(X_extend);                   % Number of the dipoles in the x-direction 
Ny=length(Y_extend);                   % Number of the dipoles in the y-direction 
Nz=length(Z_extend);                   % Number of the dipoles in the z-direction 

N=Nx*Ny*Nz;

X=reshape(x_rec,[N,1]);     % X-coordinates of the dipoles
Y=reshape(y_rec,[N,1]); % Y-coordinates of the dipoles
Z=reshape(z_rec,[N,1]);     % Z-coordinates of the dipoles

r_block=[X Y Z];         % Position of the each dipoles inside extended rectangular block
end


%将下面的复制到%%%%处可以控制平面大小相同
% distance = 1000;
% imagesize = 3200;
% scale = 1;
%     if GPU==1
%         X1=gpuArray(-round(distance/dx):round(distance/dx)); 
%         X_extend=X1*dx;
%         Y1=gpuArray(-round(0.5*scale*Max_y/dy+scale*Ly/(2*dy)):round(0.5*scale*Max_y/dy+scale*Ly/(2*dy)));
%         Y_extend=Y1*dy;
%         Z1=gpuArray(-round(0.5*scale*Max_z/dz+scale*Lz/(2*dz)):round(0.5*scale*Max_z/dz+scale*Lz/(2*dz)));
%         Z_extend=Z1*dz;
% 
%     elseif GPU==0
%         X1=-round(distance/dx):round(distance/dx); 
%         X_extend=X1*dx;
%         Y1=-round(imagesize/dy):round(imagesize/dy);
%         Y_extend=Y1*dy;
%         Z1=-round(imagesize/dy):round(imagesize/dy);
%         Z_extend=Z1*dz;
% %%%%%%%
% %控制屏幕大小相同
% % if GPU==1
% %         X1=gpuArray(-round(2*distance/(2*dx)):round(2*distance/(2*dx))); 
% %         X_extend=X1*dx;
% %         Y1=gpuArray(-round(imagesize/dy):round(imagesize/dy));
% %         Y_extend=Y1*dy;
% %         Z1=gpuArray(-round(imagesize/dy):round(imagesize/dy));
% %         Z_extend=Z1*dz;
% %         
% %         
% %     elseif GPU==0
% %         X1=-round(2*distance/(2*dx)):round(2*distance/(2*dx)); 
% % 
% %         X_extend=X1*dx;
% %         Y1=-round(0.5*scale*Max_y/dy+scale*Ly/(2*dy)):round(0.5*scale*Max_y/dy+scale*Ly/(2*dy));
% %         Y_extend=Y1*dy;
% %         Z1=-round(0.5*scale*Max_z/dz+scale*Lz/(2*dz)):round(0.5*scale*Max_z/dz+scale*Lz/(2*dz));
% %         Z_extend=Z1*dz;
% 
%     end
% 
%     [Zex,Yex]=meshgrid(Z_extend,Y_extend);
%     SIZE=size(Yex);
%     if GPU==1
%         Xex=zeros(SIZE(1,1),SIZE(1,2),'gpuArray');
%     else
%         Xex=zeros(SIZE(1,1),SIZE(1,2));
%     end
%     N_plane=SIZE(1,1)*SIZE(1,2);   % Total number of voxeles within NPs and in outside region
% 
%     %Nxex=length(X_extend);    % Number of voxels in the extended x-axis
%     %Nyex=Ny;
%     %Nzex=length(Z_extend);    % Number of voxels in the extended z-axis
%     x_plane=reshape(Xex,[N_plane,1]);   % X coordinate of the voxels in the extended x-axis
%     y_plane=reshape(Yex,[N_plane,1]);   % y coordinate of the voxels in the extended y-axis
%     z_plane=reshape(Zex,[N_plane,1]);   % z coordinate of the voxels in the extended z-axis
%     [y_rec,x_rec,z_rec]=meshgrid(Y_extend,X_extend,Z_extend);

