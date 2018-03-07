function [L] = plotVectors(vcts)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[dim,Nv] = size(vcts);
if dim ~= 2
    fprintf('Sorry, by now I only know how to deal with 2D vectors.\n')
else
    L = gobjects(1,Nv);
    for cvt = 1:Nv
        tip = vcts(:,cvt)/2;
        org = -vcts(:,cvt)/2;
        L(cvt) = plot([org(1),tip(1)],[org(2),tip(2)],'LineWidth',5);
    end
end
end

