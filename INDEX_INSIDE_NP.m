%=========================================================================%
%===== Finding index of dipoles inside NPs & ingoring other elements =====%寻找NPs内部和进入其他元素的偶极子指数
%=========================================================================%
%包含粒子的长方体，包含偶极子和空白空间两部分，这个函数的作用就是确定哪些位置是粒子——偶极子索引
function [INDEX_INSIDE]=INDEX_INSIDE_NP(GPU,X,Y,Z,N,Np_shape,Lx,Ly,Lz,Structure,d_inter,arrangement,theta)
%输入参数为GPU标志，粒子的XYZ坐标，总的粒子数，粒子的形状，包含粒子的长方体的长宽高，结构标志，双分子间距，网格标志
%========================== Monomeric structure单分子
%=========================================================================%
global orientation
if Structure=="monomeric"%单分子
    X0=0;
    Y0=0;
    Z0=0;
    
    if Np_shape=="spherical"    % NPs are sphere球形
        % Index of cubes inside NPNP内立方体的索引
        Index_in=find((sqrt((X-X0).^2+(Y-Y0).^2+(Z-Z0).^2)<=(Lx/2)));
        %找到位置到中心的距离小于半径的点，进行标记
    elseif Np_shape=="ellipsoid" % NPs are ellipsoid, head-tail orientation in z-direction
        Index_in=find((sqrt((X-X0).^2/((Lx/2)^2)+(Y-Y0).^2/((Ly/2)^2)+...
            (Z-Z0).^2/((Lz/2)^2))<=1));
    elseif Np_shape=="rod" % NPs are Rod with caps, vertically orientedNP是带帽的杆，垂直定向
        % Index of cubes inside NP   NP内立方体的索引
        Index_in=find((sqrt((X-X0).^2+(Y-Y0).^2+(abs(Z-Z0)-Lz/2+Lx/2).^2)<=(Lx/2)...
            & abs(Z-Z0)>(Lz/2-Lx/2))|((sqrt((X-X0).^2+(Y-Y0).^2)<=(Lx/2))...
            &abs(Z-Z0)<=(Lz/2-Lx/2)));
        
        
    elseif Np_shape=="rec_block" % NPs are Rectangular block, vertically oriented
        % Index of cubes inside NP
        Index_in=find((abs(X-X0)<=Lx/2 & abs(Y-Y0)<=Ly/2 & abs(Z-Z0)<=Lz/2));
    end
    %=========================================================================%
    
elseif Structure=="dimeric"
    if arrangement=="x_orient"
        X10=-(d_inter/2+Lx/2);
        Y10=0;
        Z10=0;
        X20=(d_inter/2+Lx/2);
        Y20=0;
        Z20=0;
    elseif arrangement=="y_orient"
        X10=0;
        Y10=-(d_inter/2+Ly/2);
        Z10=0;
        X20=0;
        Y20=(d_inter/2+Ly/2);
        Z20=0;
    elseif arrangement=="z_orient"
        X10=0;
        Y10=0;
        Z10=-(d_inter/2+Lz/2);
        X20=0;
        Y20=0;
        Z20=(d_inter/2+Lz/2);
    end
    
    if Np_shape=="spherical"       % NPs are sphere
        % Index of cubes inside NP
        Index_in=find((sqrt((X-X10).^2+(Y-Y10).^2+(Z-Z10).^2)<=(Lx/2))|...
            (sqrt((X-X20).^2+(Y-Y20).^2+(Z-Z20).^2)<=(Lx/2)));
        
    elseif Np_shape=="ellipsoid"   %NPs are ellipsoid and have head to tail orientation in z-direction
        % Index of cubes inside NP
        Index_in=find((sqrt((X-X10).^2/((Lx/2)^2)+(Y-Y10).^2/((Ly/2)^2)+...
            (Z-Z10).^2/((Lz/2)^2))<=1)|(sqrt((X-X20).^2/((Lx/2)^2)+...
            (Y-Y20).^2/((Ly/2)^2)+(Z-Z20).^2/((Lz/2)^2))<=1));
        
    elseif Np_shape=="rod"   %NPs are Rod with caps, vertically oriented
        % Index of cubes inside NP
        Index_in=find((sqrt((X-X10).^2+(Y-Y10).^2+(abs(Z-Z10)-Lz/2+Lx/2).^2)...
            <=(Lx/2)& abs(Z-Z10)>(Lz/2-Lx/2))|...
            ((sqrt((X-X10).^2+(Y-Y10).^2)<=(Lx/2))&abs(Z-Z10)<=(Lz/2-Lx/2))|...
            (sqrt((X-X20).^2+(Y-Y20).^2+(abs(Z-Z20)-Lz/2+Lx/2).^2)<=(Lx/2)...
            & abs(Z-Z20)>(Lz/2-Lx/2))|((sqrt((X-X20).^2+(Y-Y20).^2)<=(Lx/2))...
            &abs(Z-Z20)<=(Lz/2-Lx/2)));
       
        
        elseif Np_shape=="rec_block"     %NPs are Rectangular block, vertically oriented
        % Index of cubes inside NP
        Index_in=find((abs(X-X10)<=Lx/2 & abs(Y-Y10)<=Ly/2 & abs(Z-Z10)<=Lz/2)|...
            (abs(X-X20)<=Lx/2 & abs(Y-Y20)<=Ly/2 & abs(Z-Z20)<=Lz/2));
        
    end
end
%=========================================================================%
if GPU==1
    INDEX_INSIDE=zeros(N,1,'gpuArray'); % Defining a zero array in GPU
elseif GPU==0
    INDEX_INSIDE=zeros(N,1);            % Defining a zero array in CPU
end
% Considering dipoles inside the Np and ignoring contribution of the other terms考虑Np内部的偶极子，忽略其他项的贡献
INDEX_INSIDE(Index_in,1)=1;
end
%INDEX_INSIDE格式是[m,1]，m是偶极子数量，（x，1）中存放的值是第m个坐标位置是偶极子