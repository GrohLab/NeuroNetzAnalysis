function [feats] = getWaveformFeatures(mean_wf, fs)
%GETWAVEFORMFEATURES computes the features of the 
%   Detailed explanation goes here
tmPts = getWaveformCriticalPoints(mean_wf, fs);
[Nt, Ncl] = size(mean_wf);
tx = (0:Nt-1)/fs;
pt_dt = zeros(Ncl,1);
for ccl = 1:Ncl
    % Time difference between the peak and the trough (peak-trough time
    % diff)
    Nzc = size(tmPts{ccl,1},1);
    cwave = mean_wf(:,ccl) - mean(mean_wf(:,ccl));
    Ywf = interp1(tx,cwave,tmPts{ccl,1});
    [~, ySub] = sort(abs(Ywf));
    pt_dt(ccl) = abs(tmPts{ccl,1}(ySub(Nzc)) - tmPts{ccl,1}(ySub(Nzc-1)));
end
evals = getSpectralEntropy(mean_wf);
feats = [pt_dt, evals];
end