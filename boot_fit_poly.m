function [wmdl, mserr_fin] = boot_fit_poly(pts, n, ptsPer, it)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Npts = size(pts,1);
cpts = round(Npts * ptsPer);
mserr = zeros(it,1); cmdl = zeros(n+1, it);
pmdl = fit_poly(pts(:,1), pts(:,2), n);
[~, mserr_pop] = logCostFunction(pts, pmdl, n);
for cit = 1:it
    selSubs = randperm(Npts, cpts);
    cmdl(:, cit) = fit_poly(pts(selSubs,1), pts(selSubs,2), n);
    [~, mserr(cit)] = logCostFunction(pts, cmdl(:, cit), n);
end
mw = 1./(mserr.*sum(1./mserr)); wmdl = cmdl*mw;
[~, mserr_fin] = logCostFunction(pts, wmdl, n);
end

function [err, mserr] = logCostFunction(pts, mdl, n)
err = log(abs(pts(:,2) - pts(:,1).^(n:-1:0) * mdl) + 1);
mserr = (err'*err)/size(pts,1);
end