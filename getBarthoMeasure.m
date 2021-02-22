function [bm] = getBarthoMeasure(mean_wf, fs)
%getBarthoMeasure computes the features of the Bartho 2014 paper based on
%the waveform.
%   Detailed explanation goes here

[Nt, Ncl] = size(mean_wf);
% Window multiplication to ensure we look at the center part of the
% waveform
winMult = hann(Nt)'; 
mean_wf = mean_wf .* winMult';

tx = (0:Nt-1)'./fs - Nt/(2*fs);
tcp = getWaveformCriticalPoints(mean_wf, fs);
tcp = cellfun(@(x) x - Nt/(2*fs), tcp, 'UniformOutput', 0);
b50 = min(mean_wf) + range(mean_wf)./2;
ampCp = cell(Ncl,1);
tpd = zeros(Ncl,1); bhad = tpd;
for ccl = 1:Ncl
    ampCp(ccl) = {interp1(tx, mean_wf(:,ccl), tcp{ccl,1}, 'spline')};
    [~, troughSub] = max(abs(ampCp{ccl}));
    [~, secPeakSub] = max(abs(ampCp{ccl}(troughSub+1:end)));
    if isempty(secPeakSub)
        % We have to think of a math solution for the so called 'stable
        % region' or 'period'.
        fprintf(1, 'Cluster %d is a bitch\n', ccl)
        figure(ccl); plot(tx, mean_wf(:,ccl))
    end
    tpd(ccl) = diff(tcp{ccl,1}([troughSub, troughSub+secPeakSub]));
    if ampCp{ccl}(troughSub) > 0
        halfFlags = mean_wf(:,ccl) >= b50(ccl);
    else
        halfFlags = mean_wf(:,ccl) <= b50(ccl);
    end
    frstSub = find(halfFlags, 1, 'first');
    lstSub = find(halfFlags, 1, 'last');
    mdl = fit_poly(tx([frstSub-1, frstSub]), mean_wf([frstSub-1, frstSub],ccl),...
        1); frstXg = (b50(ccl) - mdl(2))/mdl(1);
    mdl = fit_poly(tx([lstSub-1, lstSub]), mean_wf([lstSub-1, lstSub],ccl),...
        1); lstXg = (b50(ccl) - mdl(2))/mdl(1);
    bhad(ccl) = lstXg - frstXg;
end
bm = [tpd, bhad];
end

%     if tcp{ccl, troughSub} - Ncl/(2 * fs) > 0
%         % We consider the maximum as the trough and the minimum as the peak
%         
%     end
