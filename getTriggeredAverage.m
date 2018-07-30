function [triggeredAverageSignals, signalVariation, timeAxis] =...
    getTriggeredAverage(sigStack, kIdx, timeLapse)
%GETTRIGGEREDAVERAGE Summary of this function goes here
%   Detailed explanation goes here
[~, Nt, Na] = size(sigStack);
% Time axis
timeAxis = seconds(linspace(-timeLapse(1),timeLapse(2),Nt));
if Na-sum(kIdx) >= 1
    triggeredAverageSignals = mean(sigStack(:,:,~kIdx),3);
    signalVariation = std(sigStack(:,:,~kIdx),[],3);
else
    disp('All triggered segments are set to be deleted!')
    triggeredAverageSignals = [];
    signalVariation = [];
end