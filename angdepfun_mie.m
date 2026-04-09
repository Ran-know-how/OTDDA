function [pin , taun] = angdepfun_mie(theta,n)
%ANGDEPFUN_MIE  计算角相关函数

pin = zeros(numel(n),numel(theta));
taun = zeros(numel(n),numel(theta));
mu = cos(theta);
pin(1,:) = 1;
pin(2,:) = 3*mu;
taun(1,:) = mu;
taun(2,:) = 6*mu.^2-3;
for in = 3:n(end)
    pin(in,:) = (2*in - 1)/(in-1)*mu.*pin(in-1,:)-in/(in-1).*pin(in-2,:);
    taun(in,:) = in*mu.*pin(in,:) - (in+1)*pin(in-1,:);
end
end

