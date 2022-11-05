function [bm] = getBarthoMeasure(mean_wf, fs)
%getBarthoMeasure computes the features of the Bartho 2014 paper based on
%the waveform. With some additional measures: peak-asymmetry.
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
tpd = nan(Ncl,1); bhad = tpd; asym = nan(Ncl,1);
for ccl = 1:Ncl
    ampCp(ccl) = {interp1(tx, mean_wf(:,ccl), tcp{ccl,1}, 'spline')};
    [~, troughSub] = max(abs(ampCp{ccl}));
    [~, secPeakSub] = max(abs(ampCp{ccl}(troughSub+1:end)));
    % Calculating the z-score for the found peaks
    zPeaks = (ampCp{ccl} - mean(mean_wf(:,ccl)))/std(mean_wf(:,ccl)); 
    zigFlags = abs(zPeaks) >= 1; zigPeaks = zPeaks(zigFlags);
    if sum(zigFlags) == 3
        % Looking at 'M'-shape waveforms with three critical points
        asym(ccl) = diff(zigPeaks([1,3]))/max(zigPeaks([1,3]));
    elseif sum(zigFlags) == 2
        % Looking at 'N'-shape waveforms with two critical points
        asym(ccl) = sign(diff(zigPeaks));
    end
    if isempty(secPeakSub)
        % We have to think of a math solution for the so called 'stable
        % region' or 'period'.
        fprintf(1, 'Cluster %d is a bitch\n', ccl)
        figure(ccl); plot(tx, mean_wf(:,ccl))
        continue;
    end
    tpd(ccl) = diff(tcp{ccl,1}([troughSub, troughSub+secPeakSub]));
    % if ampCp{ccl}(troughSub) > 0
    halfFlags = mean_wf(:,ccl) >= b50(ccl);
    %else
    %    halfFlags = mean_wf(:,ccl) <= b50(ccl);
    %end
    % halfFlags = mean_wf(:,ccl) >= b50(ccl);
    threshCross = diff(halfFlags);
    if sum(threshCross < 0) > 1
        fprintf(2, 'Cluster %d has wobbles\nDon''t know how to deal with this\n', ccl)
        fprintf(2, 'Skipping it...\n')
        continue
    end
    frstSub = find(threshCross, 1, 'first');
    lstSub = find(threshCross, 1, 'last');
    mdl = fit_poly(tx([frstSub, frstSub+1]), mean_wf([frstSub, frstSub+1],ccl),...
        1); frstXg = (b50(ccl) - mdl(2))/mdl(1);
    mdl = fit_poly(tx([lstSub-1, lstSub]), mean_wf([lstSub-1, lstSub],ccl),...
        1); lstXg = (b50(ccl) - mdl(2))/mdl(1);
    bhad(ccl) = lstXg - frstXg;
end
bm = [tpd, bhad, asym];
end

%     if tcp{ccl, troughSub} - Ncl/(2 * fs) > 0
%         % We consider the maximum as the trough and the minimum as the peak
%         
%     end
