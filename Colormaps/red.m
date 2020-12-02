function cmap = red(lvls)
%RED creates a colormap going from white, passing through yellow and
%arriving to red
%  
thrd = lvls/3;
sxth = thrd/2;
g = zeros(lvls,1); b = g;
% Blue channel
bSubs = (round(sxth):round(thrd))';
% Linear decay for the blue channel
bmdl = fit_poly(bSubs,linspace(1,0,numel(bSubs)),1);
b(1:bSubs(1)) = 1;
b(bSubs) = [bSubs,ones(size(bSubs))]*bmdl;
% Green channel
g(1:bSubs(end)) = 1; a = 0.01;
gSubs = (bSubs(end):lvls)';
% Exponential decay for the green channel
g(gSubs) = 1 - exp(a*(gSubs - lvls)).*linspace(0,1,numel(gSubs))';
cmap = [ones(size(g)),g,b];
cmap = cmap + min(cmap(:)); cmap = cmap./max(cmap(:));
end