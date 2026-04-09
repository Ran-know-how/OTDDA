%=========================================================================%
%===== Finding dimension of the Rectangular block that encompass NP ======%
%求包含NP的矩形块的维数
%=========================================================================%
function [Lx,Ly,Lz,theta]= Nps_parameters(r_eff,Np_shape)
%输入参量r_eff为有效尺寸，Np_shape为粒子形状

if Np_shape=="spherical"   %如果为球形   
    d_eff=2*r_eff;     % Effective diameter有效直径，半径乘2
    Lx=d_eff;
    Ly=d_eff;
    Lz=d_eff;
   %包围住粒子的尺寸就是球体的外界正方体
elseif Np_shape=="ellipsoid" %如果为椭球
    ARyx=input('\n\nEnter the ratio of y-semi axis to x-semi axis:输入y-半轴与x-半轴的比例(小数)：');
    ARzx=input('\nEnter the ratio of z-semi axis to x-semi axis:输入z半轴与x半轴的比例（小数）：');
    %以x半轴为基准，得到y半轴和z半轴
    clc
    a=((1/(ARyx*ARzx))^(1/3))*r_eff;    % Semi- minor axis in x- direction  x方向的半短轴
    b0=ARyx*a;                           % Semi- minor axis in y- direction  y方向的半短轴
    c0=ARzx*a;                           % Semi- major axis in z- direction  z方向的半短轴

    %通过有效半径的定义反推出三个半轴的长度（有效半径的定义是与粒子体积相等的球体的半径）
    Lx=2*a;                             % Length of Np in x-direction
    Ly=2*abs(b0);                             % Length of Np in y-direction
    Lz=2*abs(c0);                             % Length of Np in z-direction
    %包围住椭球的立方体长宽高就是三个方向轴的长度
elseif Np_shape=="rod"%如果是棒
    volume=4*pi/3*(r_eff^3);%棒的体积
    AR=input('\nEnter the aspect ratio, ratio of the long axis to short one:输入纵横比，长轴与短轴的比例');
    clc
    r=(volume/(pi*(2*(AR-1)+4/3)))^(1/3); % Raduis of the Rod半径
    diameter=2*r;                         % Diameter of the Rod直径
    higth=AR*diameter;                    % Higth of the Rod高度
    
    Lx=diameter;                       % Length of Np in x-direction
    Ly=diameter;                       % Length of Np in y-direction
    Lz=higth;                          % Length of Np in z-direction
    %包围住圆柱棒的立方体，高相等，圆面用正方形
elseif Np_shape=="rec_block"   %立方体           
    volume=4*pi/3*(r_eff^3);%体积
    ARyx=input('\nEnter the ratio of the y-axis to x one:输入y轴与x轴的比例');
    ARzx=input('\nEnter the ratio of the z-axis to x one:输入z轴与x轴的比例');
    clc
    a=(volume/(ARyx*ARzx))^(1/3);    % First side of rectangular block矩形块的第一侧
    b=ARyx*a;                        % Second side of rectangular block矩形块的第二侧
    c=ARzx*a;                        % Third side of rectangular block矩形块的第三侧
    %通过有效半径的尺寸反推出
    Lx=a;                            % Length of Np in x-direction
    Ly=b;                            % Length of Np in y-direction
    Lz=c;                            % Length of Np in z-direction
    %保住立方体的就是它本身
end
end