function [wmdl, inln] = boot_fit_poly(pts, n, ptsPer, it, th)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Npts = size(pts,1);
cpts = round(Npts * ptsPer);
mserr = zeros(it,1); cmdl = zeros(n+1, it);
% pmdl = fit_poly(pts(:,1), pts(:,2), n);
% [inln, mserr_pop] = meansqrerror(pts, pmdl, n);
inlnit = false(size(pts,1), it);
for cit = 1:it
    selSubs = randperm(Npts, cpts);
    cmdl(:, cit) = fit_poly(pts(selSubs,1), pts(selSubs,2), n);
    [inlnit(:,cit), mserr(cit)] = meansqrerr(pts, cmdl(:, cit), n, th);
end
consensus = mean(inlnit,2) > 0.99;
wmdl = fit_poly(pts(consensus, 1), pts(consensus, 2), n);
% mw = 1./(mserr.*sum(1./mserr)); wmdl = cmdl*mw;
inln = meansqrerr(pts, wmdl, n, th);
end

function [inln, mserr] = meansqrerr(pts, mdl, n, th)
err = pts(:,2) - ((pts(:,1).^(n:-1:0)) * mdl);
% Starting with z-score in between normal 1 std ~0.63. Later on could be a
% user defined threshold.
inln = abs(zscore(err)) <= th;
mserr = (err'*err)/size(pts,1);
end