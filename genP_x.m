function [ p_x, origSum ] = genP_x( params, x_axis)
%genP_x Generates the estimated pdf for the data from emforgmm
%   [ p_x, origSum ] = genP_x( params, x_axis)
[M,~] = size(params);
if iscolumn(x_axis)
    x_axis = x_axis';
end
p_x = zeros(1,length(x_axis));
for m = 1:M
    p_x = p_x + params(m,1)*evalgauss(x_axis,params(m,2),params(m,3));
end
origSum = sum(p_x(:));
p_x = p_x./sum(p_x(:));
end