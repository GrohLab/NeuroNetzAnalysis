function [feats] = getWaveformFeatures(mean_wf, fs)
%GETWAVEFORMFEATURES computes the features of the 
%   Detailed explanation goes here
tmPts = getWaveformCriticalPoints(mean_wf, fs);
[Nt, Ncl] = size(mean_wf);
tx = (0:Nt-1)/fs;
pt_dt = zeros(Ncl,1);
pt_r = pt_dt;
for ccl = 1:Ncl
    % Time difference between the peak and the trough (peak-trough time
    % diff)
    Nzc = size(tmPts{ccl,1},1);
    cwave = mean_wf(:,ccl) - mean(mean_wf(:,ccl));
    Ywf = interp1(tx,cwave,tmPts{ccl,1});
    [Ywf_sort, ySub] = sort(abs(Ywf));
    pt_dt(ccl) = abs(tmPts{ccl,1}(ySub(Nzc)) - tmPts{ccl,1}(ySub(Nzc-1)));
    pt_r(ccl) = Ywf_sort(2)/Ywf_sort(1);
end
[evals, mdl, gof] = getSpectralEntropy(mean_wf, [], fs);
feats = [pt_dt, evals, mdl', gof, pt_r];
end