function [ pt ] = lineIntersect( mdl1, mdl2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
xi = (mdl2(:,2) - mdl1(:,2))./(mdl1(:,1) - mdl2(:,1));
yi = mdl1(:,1).*xi + mdl1(:,2);
pt = [xi,yi];
end

