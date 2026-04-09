%在给定边长后，输出整个格点立方体中每个点的xyz坐标
function [r_block]=block(GPU,d,Max_x,Max_y,Max_z)

dx=d;
dy=d;
dz=d;
    
    if GPU==1
        X1=gpuArray(-round(Max_x/(2*dx)):round(Max_x/(2*dx)));
        X_extend=X1*dx;
        Y1=gpuArray(-round(Max_y/(2*dy)):round(Max_y/(2*dy)));
        Y_extend=Y1*dy;
        Z1=gpuArray(-round(Max_z/(2*dz)):round(Max_z/(2*dz)));
        Z_extend=Z1*dz;
        
    elseif GPU==0
        X1=-round(Max_x/(2*dx)):round(Max_x/(2*dx));
        X_extend=X1*dx;
        Y1=-round(Max_y/(2*dy)):round(Max_y/(2*dy));
        Y_extend=Y1*dy;
        Z1=-round(Max_z/(2*dz)):round(Max_z/(2*dz));
        Z_extend=Z1*dz;
    end
    
%     [x_rec,y_rec,z_rec]=meshgrid(X_extend,Y_extend,Z_extend);
    [y_rec,x_rec,z_rec]=meshgrid(Y_extend,X_extend,Z_extend);

Nx=length(X_extend);                   % Number of the dipoles in the x-direction 
Ny=length(Y_extend);                   % Number of the dipoles in the y-direction 
Nz=length(Z_extend);                   % Number of the dipoles in the z-direction 

N=Nx*Ny*Nz;

X=reshape(x_rec,[N,1]);     % X-coordinates of the dipoles每个偶极子在扩展矩形块内的x坐标
Y=reshape(y_rec,[N,1]); % Y-coordinates of the dipoles
Z=reshape(z_rec,[N,1]);     % Z-coordinates of the dipoles

r_block=[X Y Z];
end