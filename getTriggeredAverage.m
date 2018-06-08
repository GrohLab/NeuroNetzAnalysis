function [Wmean, Wstd] = getTriggeredAverage(sigStack, kIdx)
%GETTRIGGEREDAVERAGE Summary of this function goes here
%   Detailed explanation goes here
Na = size(sigStack, 2);
if Na-sum(~kIdx) >= 1
    Wmean = mean(sigStack(:,~kIdx),2);
    Wstd = std(sigStack(:,~kIdx),[],2);
else
    disp('All triggered segments are set to be deleted!')
    Wmean = [];
    Wstd = [];
end

