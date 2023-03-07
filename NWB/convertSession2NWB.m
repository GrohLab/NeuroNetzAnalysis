function [nwbObj] = convertSession2NWB(sessionPath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

addRequired(p, "sessionPath", @(x) exist(x, "dir"));

parse(p, sessionPath);
sessionPath = p.Results.sessionPath;
%%
expandName = @(x) fullfile(x.folder, x.name);
chanMapDir = "Z:\Emilio";
chanMapPaths = dir(fullfile(chanMapDir, "Cambridge_NeuroTech_*.mat"));
chanMapPaths(arrayfun(@(x) contains(x.name, "old"), chanMapPaths)) = [];
chanMapTypes = extractBetween(arrayfun(@(x) string(x.name), ...
    chanMapPaths), "Cambridge_NeuroTech_",".mat");

[~, ephSes] = fileparts(sessionPath);
chanMapUsed = extractAfter(ephSes,"ephys_");
cmFlag = contains(chanMapTypes, chanMapUsed);
% Metadata for the electrode object.
metaData = {'BrainStructure', 'Superior colliculus',...
            'Coordinates', [1500, -3600, -3200], ...
            'Filtering', 'anti-aliasing'};
% Initialise NwbFile object
nwbObj = getSessionInfo4NWB(sessionPath, ...
    'Comment', 'Awake head-fixed; puffed from behind');
% Add silicon probe information
nwbObj = CNTP2NWB(nwbObj, expandName(chanMapPaths(cmFlag)), metaData{:});
% Subject info
subject = types.core.Subject(...
    'subject_id', extractBefore(nwbObj.identifier,"_"), ...
    'date_of_birth', datetime(), ...
    'species', 'mus musculus', ...
    'sex', 'male', ...
    'description', 'IBF');
nwbObj.general_subject = subject;

nwbObj = getTrialInfo4NWB(nwbObj, sessionPath, ...
    'ConditionSelection', 'Control Puff');


end