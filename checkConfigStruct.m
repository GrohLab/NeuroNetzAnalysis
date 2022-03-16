function [iOk, fnH] = checkConfigStruct(cStruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
checkCS = @(x) isstruct(x) & all(contains(fieldnames(x), ...
    {'Experiment', 'Viewing_window_s', 'Response_window_s', 'BinSize_s', ...
    'Trigger', 'ConsideredConditions'}));
iOk = checkCS(cStruct); fnH = checkCS;
end