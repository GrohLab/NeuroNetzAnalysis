function [corrStackCells] = neuroCorr(spks, tmReach, binSz, fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% getStacks(spT,alignP,ONOFF,timeSpan,fs,fsLFP,consEvents,...)
corrStackCells = cell(length(spks),1); Ncl = length(spks);
zeroCorr = ceil(tmReach * fs) + 1;
for ccl = 1:size(spks, 1)
    corrStackCells(ccl) = {squeeze(sum(getStacks(false, spks{ccl}, 'on',...
        [-tmReach, tmReach], fs, fs, spks(ccl+1:Ncl)),3))};
    if isempty(corrStackCells{ccl})
        Ne = size(ccl+1:Ncl,2); Nt = 2*ceil(tmReach*fs) + 1;
        fprintf(1, 'Decreasing the number of spikes for cluster %d\n',ccl)
        m = memory; Nspks = m.MaxPossibleArrayBytes / ((2+Ne) * Nt);
        Nspks_old = size(spks{ccl},1);
        fprintf(1, 'From %d to %d\n', Nspks_old, Nspks)
        randSubs = sort(randperm(Nspks_old, round(Nspks))); 
        spks{ccl} = spks{ccl}(randSubs);
        corrStackCells(ccl) = {squeeze(sum(getStacks(false, spks{ccl}, 'on',...
            [-tmReach, tmReach], fs, fs, spks(ccl+1:Ncl)),3))};
    end
    corrStackCells{ccl}(2,:) = []; corrStackCells{ccl}(1,zeroCorr) = 0;
end
end
% corrStackCells = cellfun(@(x) squeeze(sum(getStacks(false, x, 'on',...
%     [-tmReach, tmReach], fs, fs, spks),3)), spks, 'UniformOutput', 0);