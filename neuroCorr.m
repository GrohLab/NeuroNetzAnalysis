function [corrStackCells] = neuroCorr(spks, tmReach, binSz, fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% getStacks(spT,alignP,ONOFF,timeSpan,fs,fsLFP,consEvents,...)
corrStackCells = cell(length(spks),1); Ncl = length(spks);
zeroCorr = ceil(tmReach * fs) + 1;
for ccl = 1:size(spks, 1)
    corrStackCells(ccl) = {squeeze(sum(getStacks(false, spks{ccl}, 'on',...
        [-tmReach, tmReach], fs, fs, spks(ccl+1:Ncl)),3))};
    corrStackCells{ccl}(2,:) = []; corrStackCells{ccl}(:,zeroCorr) = 0;
end
end
% corrStackCells = cellfun(@(x) squeeze(sum(getStacks(false, x, 'on',...
%     [-tmReach, tmReach], fs, fs, spks),3)), spks, 'UniformOutput', 0);