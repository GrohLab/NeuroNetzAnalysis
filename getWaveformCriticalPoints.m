function timePts = getWaveformCriticalPoints(avWaves, fs)
% GETWAVEFORMCRITICALPOINTS returns the time points for which the first and
% second derivatives of the signal equals zero for the minima and/or maxima
% and the inflection points.

dw = diff(avWaves, 1, 1);
ddw = diff(avWaves, 2, 1);
dt = 1/fs;
[Nt, Ncl] = size(avWaves);
tx = (0:Nt-1)/fs;
timePts = cell(Ncl,2);
for ccl = 1:Ncl
    zcIdx = diff(sign(dw(:,ccl))) ~= 0;
    fpIdx = diff(sign(ddw(:,ccl))) ~= 0;
    zcSubs = find(zcIdx);
    fpSubs = find(fpIdx);
    xz = zeros(numel(zcSubs),1);
    for czc = 1:numel(zcSubs)
        lnIdx = zcSubs(czc):zcSubs(czc)+1;
        mdl = fit_poly(tx(lnIdx), dw(lnIdx, ccl), 1);
        xz(czc) = -mdl(2)/mdl(1) + dt/2;
    end
    
    fp = zeros(numel(fpSubs),1);
    for cfp = 1:numel(fpSubs)
        lnIdx = fpSubs(cfp):fpSubs(cfp)+1;
        mdl = fit_poly(tx(lnIdx), ddw(lnIdx, ccl), 1);
        fp(cfp) = -mdl(2)/mdl(1) + dt;
    end
%     lmt = abs(median(devs(ccl,:)));
%     xzImp = abs(interpolateTimeSeries(tx, devs(ccl,:), xz)) > lmt;
%     fpImp = abs(interpolateTimeSeries(tx, devs(ccl,:), fp)) > lmt;
    timePts(ccl,:) = {xz, fp};
end


end