function [nwbObj, trial_timeseries, configStructure] = ...
    getTrialInfo4NWB(nwbObj, sessionPath, varargin)
%getTrialInfo4NWB Summary of this function goes here
%   [nwbObj, trial_timeseries, configStructure] = getTrialInfo4NWB(nwbObj, sessionPath)
%   ... = getTrialInfo4NWB(..., 'Name','Value')
%       INPUTS:
%           nwbObj - NwbFile object
%           sessionPath - string or char indicating the folder with all
%           ephys files for a session on the roller.
%           Name-Value:
%               ConditionSelection - string, char or cellstr variable
%               indicating the condition names to consider in the session
%       OUTPUTS:
%           nwbObj - NwbFile object
%           trial_timeseries - cell array with necessary information per
%           trial for NWB framework. Will be used in |assignSpkTms2NWB|
%           function.
%           configStructure - structure with analysis metadata for the
%           session.
%Emilio Isaias-Camacho @GrohLab 2023
%% Validate inputs
p = inputParser();

istxt = @(x) isstring(x) | ischar(x) | iscellstr(x);
validateCondSel = @(x) isPositiveIntegerValuedNumeric(x) | istxt(x);

addRequired(p, "nwbObj", @(x) isa(x, "NwbFile"))
addRequired(p, "sessionPath", @(x) exist(x, 'dir'))
addParameter(p, "ConditionSelection", "all", validateCondSel)

parse(p, nwbObj, sessionPath, varargin{:})

nwbObj = p.Results.nwbObj;
sessionPath = p.Results.sessionPath;
condSel = p.Results.ConditionSelection;

%% Auxiliary functions and variables
fnOpts = {'UniformOutput', false};
pathHere = @(x) fullfile(sessionPath, x);
look4this = @(x) dir(pathHere(x));
expandName = @(x) fullfile(x.folder, x.name);
% Validating for files in the ephys directory
cFiles = look4this("*analysis.mat");
fFiles = look4this("*_sampling_frequency.mat");
pFiles = look4this("*_ParameterConfiguration.mat");
if isempty(cFiles) || isempty(fFiles) || isempty(pFiles)
    fprintf(1, "Warning! Necessary files missing in given directory!\n")
    fprintf(1, "Cannot continue!\n")
    return
end
cMF = matfile(expandName(cFiles));
Conditions = cMF.Conditions; Triggers = cMF.Triggers; clearvars("cMF");
condNames = arrayfun(@(x) string(x.name), Conditions);
allCond = contains(condNames, "all", "IgnoreCase", true);
if istxt(condSel)
    if strcmp(condSel, "all")
        % All conditions into the NWB file
        condSel = setdiff(1:numel(Conditions),find(allCond));
    else
        % Selected conditions by their name or it might be a cell array
        % with char arrays or strings
        condSel = ismember(condNames, condSel);
    end
elseif any(condSel-sum(allCond) > sum(~allCond))
    % in development.
end

%% Loading and preparing data
% Triggers per selected condition
trigPerCond = arrayfun(@(x) x.Triggers(:,1), Conditions(condSel), ...
    fnOpts{:});
% Number of trials per condition
Na = cellfun(@(c) size(c,1), trigPerCond);
% Trial types: Condition names
consCondNames = arrayfun(@(c) string(c.name), Conditions(condSel));
load(expandName(fFiles), "fs")
load(expandName(pFiles), "configStructure")

csFlag = arrayfun(@(x) ismember(x.ConsideredConditions, ...
    consCondNames), configStructure, fnOpts{:});
csFlag = find(cellfun(@any, csFlag),1,"first");
configStructure = configStructure(csFlag); clearvars("csFlag");
offst = configStructure.Viewing_window_s(1);
configStructure.Sampling_rate_hz = fs;
tsNames = string(fieldnames(Triggers));
% Trial duration
trial_duration = diff(configStructure.Viewing_window_s);
trial_samples = ceil(trial_duration*fs);
% Trial start in samples
trial_start_sample = cellfun(@(c) c + ...
    round(configStructure.Viewing_window_s(1)*fs), trigPerCond, fnOpts{:}); 
ccfn = 1;
N_signals = numel(tsNames); N_trial_types = numel(consCondNames);
chopped_timeseries = cell(N_signals, N_trial_types);
% Chopping signals around the triggers
for cfn = tsNames(:)'
    sgn = Triggers.(cfn);
    chopped_timeseries_aux = cellfun(@(c) arrayfun(@(t) ...
        sgn(t:t+trial_samples-1), c, fnOpts{:}), ...
        trial_start_sample, fnOpts{:});
    chopped_timeseries_aux = cellfun(@(st) cat(1, st{:}), ...
        chopped_timeseries_aux, fnOpts{:});
    chopped_timeseries(ccfn,:) = chopped_timeseries_aux;
    ccfn = ccfn + 1;
end
% Freeing memory from auxiliary variables
clearvars("sgn", "chopped_timeseries_aux")
% Trial start in time
trial_start_time_per_cond = cellfun(@(t) (t./fs) + offst, ...
    trigPerCond, fnOpts{:});
timestamps = cellfun(@(c) arrayfun(@(t) t:t+trial_samples-1, c, ...
    fnOpts{:}), trial_start_sample, fnOpts{:});
timestamps = cellfun(@(c) cat(2, c{:})./fs, timestamps, fnOpts{:});
timestamps_all = cat(2, timestamps{:});
trial_type = arrayfun(@(tt, tm) repmat(string(tt), tm, 1), ...
    consCondNames, Na, fnOpts{:}); trial_type = cat(1, trial_type{:});
%% Variable storing continuous signals from the trials.
% Trial timeseries

timeseries_aux = cellfun(@(c) reshape(transpose(c), numel(c), []), ...
        chopped_timeseries, fnOpts{:});
timeseries_aux = arrayfun(@(x) single(cat(1, timeseries_aux{x,:})), ...
    1:N_signals, fnOpts{:});

ts_ref = arrayfun(@(x) types.untyped.ObjectView([...
    '/stimulus/presentation/' char(x)]), tsNames);
for cfn = tsNames(:)'
    ccfn = ismember(tsNames, cfn);
    if strcmpi(cfn,'laser')
        ts_aux = types.core.OptogeneticSeries();
        ts_aux.site = types.core.OptogeneticStimulusSite(...
                'description', 'Laser fiber attached to the electrode',...
                'device', types.untyped.SoftLink(nwbObj.general_devices.get('Laser')),...
                'excitation_lambda', 488,...
                'location', 'Electrode tip, in superior colliculus SCig');
    else
        ts_aux = types.core.TimeSeries();
    end
    ts_aux.data = timeseries_aux{ccfn};
    ts_aux.data_unit = 'ADC_16bits';
    ts_aux.description = [char(cfn) ' stimulation'];
    ts_aux.timestamps = timestamps_all(:);
    
    nwbObj.stimulus_presentation.set(cfn, ts_aux);
end

trial_timeseries = arrayfun(@(x) {ts_ref(1) int64(x) int64(trial_samples); ...
               ts_ref(2) int64(x) int64(trial_samples)}, ...
               cat(1, trial_start_sample{:}), fnOpts{:});
%% Trial epoch - from tutorial
% Though TimeIntervals is a subclass of the DynamicTable type, we opt for
% populating the Dynamic Table data by column instead of using `addRow`
% here because of how the data is formatted. DynamicTable is flexible
% enough to accomodate both styles of data conversion.
trials_epoch = types.core.TimeIntervals(...
    'start_time', types.hdmf_common.VectorData('description', ...
        'Start time of trial, in seconds.'),...
    'stop_time', types.hdmf_common.VectorData('description', ...
        'Stop time of trial, in seconds.'),...
    'colnames', [cellstr(consCondNames);...
        {'start_time'; 'stop_time'; 'timeseries'}],...
    'description', 'trial data and properties', ...
    'id', types.hdmf_common.ElementIdentifiers(),...
    'timeseries', types.core.TimeSeriesReferenceVectorData('description', 'Index into timeseries Data'),...
    'timeseries_index', types.hdmf_common.VectorIndex(...
        'description', 'Index into Timeseries VectorData',...
        'target', types.untyped.ObjectView('/intervals/trials/timeseries')));

trials_epoch.start_time.data = cat(1, trial_start_time_per_cond{:});
trials_epoch.stop_time.data = cat(1, trial_start_time_per_cond{:}) + ...
    trial_duration;
trials_epoch.id.data = (1:sum(Na))';
arrayfun(@(cn) trials_epoch.vectordata.set(cn, ...
    types.hdmf_common.VectorData('data', contains(trial_type, cn),...
    'description', char(cn))), consCondNames);

nwbObj.intervals_trials = trials_epoch;
end