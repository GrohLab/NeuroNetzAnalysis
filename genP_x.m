function [ p_x ] = genP_x( params, x_axis)
%genP_x Generates the estimated pdf for the data from emforgmm
[M,~] = size(params);
p_x = zeros(1,length(x_axis));
for m = 1:M
    p_x = p_x + params(m,1)*evalgauss(x_axis,params(m,2),params(m,3));
end
p_x = p_x./sum(p_x(:));
end