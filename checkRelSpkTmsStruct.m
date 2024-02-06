function [iOk, fnH] = checkRelSpkTmsStruct(rstStruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fnCell = {'name', 'SpikeTimes'};
checkRST = @(x) isstruct(x) & (all(contains(fieldnames(x), ...
    fnCell(1:end-1))) | all(contains(fieldnames(x), fnCell)));
iOk = checkRST(rstStruct); fnH = checkRST;
end