function [iOk, fnH] = checkConfigStruct(cStruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fnCell = {'Experiment', 'Viewing_window_s', 'Response_window_s', ...
    'BinSize_s', 'Trigger', 'ConsideredConditions', ...
    'Spontaneous_window_s'};
checkCS = @(x) isstruct(x) & (all(contains(fieldnames(x), ...
    fnCell(1:end-1))) | all(contains(fieldnames(x), fnCell)));
iOk = checkCS(cStruct); fnH = checkCS;
end