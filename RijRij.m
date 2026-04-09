
%====== Calculating RijRij-I3 and 3RijRij-I3 in interaction matrix A
%=====%计算相互作用矩阵
%=========================================================================%
function [rjkrjk1_I,rjkrjk2_I,rjkrjk3_I,rjkrjk4_I,rjkrjk5_I,rjkrjk6_I,rjkrjk31_I,...
    rjkrjk32_I,rjkrjk33_I,rjkrjk34_I,rjkrjk35_I,rjkrjk36_I,RJK]=RijRij(r_block)
%输入参数：r_block，形状为[N,3],r_block(N,:)可以提取出第N个空格的xyz坐标
rkj1=r_block(1,1)-r_block(:,1);%第i个粒子到第j个粒子的x坐标差
rkj2=r_block(1,2)-r_block(:,2);%第i个粒子到第j个粒子的y坐标差
rkj3=r_block(1,3)-r_block(:,3);%第i个粒子到第j个粒子的z坐标差


rk_to_rj=[rkj1  rkj2  rkj3];%[N,3]矩阵，（N,;）可以得到第一个粒子与其余所有粒子的xyz坐标差
rk_to_rj(1,:)=1e-9;                % in order to scape from NAN Error以便从NAN错误中查找，将自己与自己的距离改为1，防止除0错误
RJK=sqrt(rk_to_rj(:,1).^2+rk_to_rj(:,2).^2+rk_to_rj(:,3).^2);%得到[N,1]的数列，表示第一个粒子与其他所有粒子的距离
rjkrjk=[rkj1./RJK  rkj2./RJK  rkj3./RJK];%[N,3]阵列，坐标差除以距离，物理意义是方向余弦，单位矢量


rjkrjk1_I=rjkrjk(:,1).*rjkrjk(:,1)-1;%计算rjkrjk第一列的平方减1
rjkrjk2_I=rjkrjk(:,1).*rjkrjk(:,2);%计算rjkrjk第一列与第二列的乘积
rjkrjk3_I=rjkrjk(:,1).*rjkrjk(:,3);%计算rjkrjk第一列与第三列的乘积
rjkrjk4_I=rjkrjk(:,2).*rjkrjk(:,2)-1;%计算rjkrjk第二列的平方减1
rjkrjk5_I=rjkrjk(:,2).*rjkrjk(:,3);%计算rjkrjk第二列与第三列的乘积
rjkrjk6_I=rjkrjk(:,3).*rjkrjk(:,3)-1;%计算rjkrjk第三列的平方减1
%上面6行计算了一个上三角矩阵减去单位矩阵之后每个元素的值
%相互作用矩阵A是一个上三角矩阵，但是这里算的显然不是A，而是A的前身
%为什么是rjk的平方呢？——并在一起运算是相乘吗？

rjkrjk31_I=3*rjkrjk(:,1).*rjkrjk(:,1)-1;%计算3倍rjkrjk第一列的平方减1
rjkrjk32_I=3*rjkrjk(:,1).*rjkrjk(:,2);%计算3倍rjkrjk第一列与第二列的乘积
rjkrjk33_I=3*rjkrjk(:,1).*rjkrjk(:,3);%计算3倍rjkrjk第一列与第三列的乘积
rjkrjk34_I=3*rjkrjk(:,2).*rjkrjk(:,2)-1;%计算3倍rjkrjk第二列的平方减1
rjkrjk35_I=3*rjkrjk(:,2).*rjkrjk(:,3);%计算3倍rjkrjk第二列与第三列的乘积
rjkrjk36_I=3*rjkrjk(:,3).*rjkrjk(:,3)-1;%计算3倍rjkrjk第三列的平方减1
%这里也是两个方向矢量并在一起乘3减去单位矢量，详细见A的公式《离散偶极子近似公式推导》

rjkrjk1_I(1,1)=0;%将第一个元素设置为0，防止除0错误，以下全是
rjkrjk2_I(1,1)=0;
rjkrjk3_I(1,1)=0;
rjkrjk4_I(1,1)=0;
rjkrjk5_I(1,1)=0;
rjkrjk6_I(1,1)=0;

rjkrjk31_I(1,1)=0;
rjkrjk32_I(1,1)=0;
rjkrjk33_I(1,1)=0;
rjkrjk34_I(1,1)=0;
rjkrjk35_I(1,1)=0;
rjkrjk36_I(1,1)=0;

end
%=========================================================================%