function [nwbObj, trial_timeseries] = ...
    getTrialInfo4NWB(nwbObj, sessionPath, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Validate inputs
p = inputParser();
validateCondSel = @(x) isPositiveIntegerValuedNumeric(x) | (isstring(x) | ...
    ischar(x));

addRequired(p, "nwbObj", @(x) isa(x, "NwbFile"))
addRequired(p, "sessionPath", @(x) exist(x, 'dir'))
addParameter(p, "ConditionSelect", "all", validateCondSel)

parse(p, nwbObj, sessionPath)

nwbObj = p.Results.nwbObj;
sessionPath = p.Results.sessionPath;

%% Auxiliary functions and variables
pathHere = @(x) fullfile(sessionPath, x);
look4this = @(x) dir(pathHere(x));
expandName = @(x) fullfile(x.folder, x.name);
% Validating for files in the ephys directory
cFiles = look4this("*analysis.mat");
if isempty(cFiles)
    fprintf(1, "Warning! No analysis file found in given directory!\n")
    fprintf(1, "Cannot continue!\n")
    return
end
cMF = matfile(expandName(cFiles));
Conditions = cMF.Conditions; Triggers = cMF.Triggers;


% Variable storing continuous signals from the trials.
trial_timeseries = cell(size(data.trialIds)); 

%%
% NWB comes with default support for trial-based data.  These must be *TimeIntervals* that
% are placed in the |intervals| property.  Note that |trials| is a special
% keyword that is required for PyNWB compatibility.

ephus = data.timeSeriesArrayHash.value{1};
ephusUnit = data.timeUnitNames{data.timeUnitIds(ephus.timeUnit)};

% lick direction and timestamps trace
tsIdx = strcmp(ephus.idStr, 'lick_trace');
bts = types.core.BehavioralTimeSeries();

bts.timeseries.set('lick_trace_ts', ...
    types.core.TimeSeries(...
    'data', ephus.valueMatrix(:,tsIdx),...
    'data_unit', ephusUnit,...
    'description', ephus.idStrDetailed{tsIdx}, ...
    'timestamps', ephus.time, ...
    'timestamps_unit', ephusUnit));
nwb.acquisition.set('lick_trace', bts);
bts_ref = types.untyped.ObjectView('/acquisition/lick_trace/lick_trace_ts');

% acousto-optic modulator input trace
tsIdx = strcmp(ephus.idStr, 'aom_input_trace');
ts = types.core.TimeSeries(...
    'data', ephus.valueMatrix(:,tsIdx), ...
    'data_unit', 'Volts', ...
    'description', ephus.idStrDetailed{tsIdx}, ...
    'timestamps', ephus.time, ...
    'timestamps_unit', ephusUnit);
nwb.stimulus_presentation.set('aom_input_trace', ts);
ts_ref = types.untyped.ObjectView('/stimulus/presentation/aom_input_trace');

% laser power
tsIdx = strcmp(ephus.idStr, 'laser_power');
ots = types.core.OptogeneticSeries(...
    'data', ephus.valueMatrix(:,tsIdx), ...
    'data_unit', 'mW', ...
    'description', ephus.idStrDetailed{tsIdx}, ...
    'timestamps', ephus.time, ...
    'timestamps_unit', ephusUnit, ...
    'site', types.untyped.SoftLink('/general/optogenetics/photostim'));
nwb.stimulus_presentation.set('laser_power', ots);
ots_ref = types.untyped.ObjectView('/stimulus/presentation/laser_power');

% append trials timeseries references in order
[ephus_trials, ~, trials_to_data] = unique(ephus.trial);
for i=1:length(ephus_trials)
    i_loc = i == trials_to_data;
    t_start = find(i_loc, 1);
    t_count = sum(i_loc);
    trial = ephus_trials(i);
    
    trial_timeseries{trial}(end+(1:3), :) = {...
        bts_ref int64(t_start) int64(t_count);...
        ts_ref  int64(t_start) int64(t_count);...
        ots_ref int64(t_start) int64(t_count)};
end
% Though TimeIntervals is a subclass of the DynamicTable type, we opt for
% populating the Dynamic Table data by column instead of using `addRow`
% here because of how the data is formatted. DynamicTable is flexible
% enough to accomodate both styles of data conversion.
trials_epoch = types.core.TimeIntervals(...
    'start_time', types.hdmf_common.VectorData('description', ...
        'Start time of epoch, in seconds.'),...
    'colnames', [data.trialTypeStr; data.trialPropertiesHash.keyNames .';...
        {'start_time'; 'stop_time'; 'acquisition'; 'timeseries'}],...
    'description', 'trial data and properties', ...
    'id', types.hdmf_common.ElementIdentifiers(),...
    'timeseries', types.core.TimeSeriesReferenceVectorData('description', 'Index into timeseries Data'),...
    'timeseries_index', types.hdmf_common.VectorIndex(...
    'description', 'Index into Timeseries VectorData',...
    'target', types.untyped.ObjectView('/intervals/trials/timeseries')));

trials_epoch.start_time.data = data.trialStartTimes;
trials_epoch.id.data = data.trialIds;

for i=1:length(data.trialTypeStr)
    trials_epoch.vectordata.set(data.trialTypeStr{i}, ...
        types.hdmf_common.VectorData('data', data.trialTypeMat(i,:),...
            'description', data.trialTypeStr{i}));
end

for i=1:length(data.trialPropertiesHash.keyNames)
    descr = data.trialPropertiesHash.descr{i};
    if iscellstr(descr)
        descr = strjoin(descr, newline);
    end
    trials_epoch.vectordata.set(data.trialPropertiesHash.keyNames{i}, ...
        types.hdmf_common.VectorData(...
        'data', data.trialPropertiesHash.value{i}, ...
        'description', descr));
end
end