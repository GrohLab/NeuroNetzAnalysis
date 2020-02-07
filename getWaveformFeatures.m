function [feats] = getWaveformFeatures(mean_wf, fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tmPts = getWaveformCriticalPoints(mean_wf, fs);
[Nt, Ncl] = size(mean_wf);
tx = (0:Nt-1)/fs;
pt_dt = zeros(Ncl,1);
for ccl = 1:Ncl
    % Time difference between the peak and the trough (peak-trough time
    % diff)
    Nzc = size(tmPts{ccl,1},1);
    Ywf = interp1(tx,mean_wf(:,ccl),tmPts{ccl,1});
    [~, ySub] = sort(Ywf);
    pt_dt(ccl) = abs(tmPts{ccl,1}(ySub(1)) - tmPts{ccl,1}(ySub(Nzc)));
end
evals = getSpectralEntropy(mean_wf);
feats = [pt_dt, evals];
end

