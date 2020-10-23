function [] = neuroCorr(spks, tmReach, binSz, fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% getStacks(spT,alignP,ONOFF,timeSpan,fs,fsLFP,consEvents,...)
maxN = max(cellfun(@max, spks)) + round(tmReach * fs);
corrStackCells = cellfun(@(x) squeeze(sum(getStacks(false(1,maxN), x, 'on',...
    [-tmReach, tmReach], fs, fs, spks),3)), spks, 'UniformOutput', 0);
end