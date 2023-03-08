function [nwbObj] = convertSession2NWB(sessionPath, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

addRequired(p, "sessionPath", @(x) exist(x, "dir"));
addParameters 
parse(p, sessionPath, varargin{:});
sessionPath = p.Results.sessionPath;
%TODO: Implement the session specific data.
%%
expandName = @(x) fullfile(x.folder, x.name);
outDir = "Z:\SC Anatomy paper data\Roller";
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
% Subject info
subject = types.core.Subject(...
    'subject_id', extractBefore(nwbObj.identifier,'_'), ...
    'date_of_birth', datetime(2021,05,20), ...
    'species', 'mus musculus', ...
    'sex', 'male', ...
    'description', 'IBF');
nwbObj.general_subject = subject;
% Add silicon probe information
nwbObj = CNTP2NWB(nwbObj, expandName(chanMapPaths(cmFlag)), metaData{:});
% Add trial information: type, start time, duration
condSel = {'ConditionSelection', 'Control Puff'};
[nwbObj, trial_timeseries, configStructure] = getTrialInfo4NWB(nwbObj, ...
    sessionPath, condSel{:});
% Flatten spike times for NWB. From cell array to 'jagged' or 'staggered'
% arrays
rstNWB = prepareRelSpkTms4NWB(sessionPath, configStructure, condSel{:});
% Add spike times per unit per trial in NWB object
[nwbObj, trial_timeseries] = assignSpkTms2NWB(nwbObj, rstNWB, ...
    trial_timeseries, configStructure);

nwbObj = exportNWB(nwbObj, trial_timeseries, 'OutDir', outDir);
end