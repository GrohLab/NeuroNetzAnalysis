function [nwbObj] = convertSession2NWB(sessionPath, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

istxt = @(x) isstring(x) | ischar(x);
% Coordinates are three dimensional; (+x = posterior, +y = inferior, 
% +z = subject s right).
validateCoords = @(x) isnumeric(x) & numel(x)==3 & isvector(x);

addRequired(p, "sessionPath", @(x) istxt(x) & exist(x, "dir"));
addParameter(p, "Birthdate", datetime, @isdatetime)
addParameter(p, "Coordinates", zeros(3,1), validateCoords)
addParameter(p, "OutDir", "Z:\SC Anatomy paper data", istxt);
addParameter(p, "ConditionSelection", "all", @(x) istxt(x) | iscellstr(x)) %#ok<ISCLSTR> 
addParameter(p, "MouseSource", 'IBF', istxt)
addParameter(p, "MouseSex", 'male', istxt)

parse(p, sessionPath, varargin{:});

sessionPath = p.Results.sessionPath;
birthdate = p.Results.Birthdate;
coords = p.Results.Coordinates;
outDir = p.Results.OutDir;
consCondName = p.Results.ConditionSelection;
mouseSource = p.Results.MouseSource;
mouseSex = p.Results.MouseSex;
%TODO: Implement the session specific data.
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
            'Coordinates', coords, ...
            'Filtering', 'anti-aliasing'};
% Initialise NwbFile object
nwbObj = getSessionInfo4NWB(sessionPath, ...
    'Comment', 'Awake head-fixed; contra-lateral whiskers puffed from behind');
% Subject info
subject = types.core.Subject(...
    'subject_id', extractBefore(nwbObj.identifier,'_'), ...
    'date_of_birth', birthdate, ...
    'species', 'mus musculus', ...
    'sex', mouseSex, ...
    'description', mouseSource);
nwbObj.general_subject = subject;
% Add silicon probe information
nwbObj = CNTP2NWB(nwbObj, expandName(chanMapPaths(cmFlag)), metaData{:});
% Add trial information: type, start time, duration
condSel = {'ConditionSelection', consCondName};
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