function [nwbObj, trial_timeseries] = assignSpkTms2NWB(nwbObj, rstNWB, ... 
    trial_timeseries, configStructure)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Input parsing
p = inputParser;

validateStructure = @(s, f) isstruct(s) & all(isfield(s, f));

rst_fields = {'TrialType', 'TrialID', 'UnitID', 'SpikeTimes', 'NumTrials',...
    'NumUnits'};
validateRST = @(x) validateStructure(x, rst_fields);

cs_fields = {'Experiment', 'Viewing_window_s', 'Response_window_s', ...
    'Spontaneous_window_s', 'BinSize_s', 'Trigger', ...
    'ConsideredConditions', 'Sampling_rate_hz'};
validateCS = @(x) validateStructure(x, cs_fields);

addRequired(p, "nwbObj", @(x) isa(x, "NwbFile"))
addRequired(p, "rstNWB", validateRST)
addRequired(p, "trial_timeseries", @iscell)
addRequired(p, "configStructure", validateCS)

parse(p, nwbObj, rstNWB, trial_timeseries, configStructure)

nwbObj = p.Results.nwbObj;
rstNWB = p.Results.rstNWB;
trial_timeseries = p.Results.trial_timeseries;
configStructure = p.Results.configStructure;
%% Auxiliary functions and variables
sessionPath = fileparts(configStructure.Experiment);
fnOpts = {'UniformOutput', false};
% expandName = @(x) fullfile(x.folder,x.name);
pathHere = @(x) fullfile(sessionPath, x);
% look4this = @(x) dir(pathHere(x));
%% Loading and/or getting waveform data
fs = configStructure.Sampling_rate_hz;
% wFiles = look4this("*_waveforms.mat");
% if isempty(wFiles)
%     clWaveforms = getClusterWaveform(rstNWB.UnitID, sessionPath);
% else
%     load(expandName(wFiles), 'clWaveforms')
%     if size(clWaveforms, 1) ~= size(rstNWB.UnitID,1)
%         clWaveforms = getClusterWaveform(rstNWB.UnitID, sessionPath);
%     end
% end
% [~, numOrd] = ismember(rstNWB.UnitID, clWaveforms(:,1));
% clWaveforms = clWaveforms(numOrd, :); clearvars("numOrd","wFiles")
%% Preparing variables per unit
% Cluster information
clInfo = getClusterInfo(pathHere("cluster_info.tsv"));
clInfo = clInfo(rstNWB.UnitID,:);
% Unit IDs
ids = uint8(str2double(rstNWB.UnitID));
% Trial ID per spike
tIDs = cat(2, rstNWB.TrialID{:}); tIDs = arrayfun(@(r) ...
    cat(2, tIDs{r,:}), (1:size(tIDs,1))', fnOpts{:});
% Trigger relative spikes
spk_events = cat(2, rstNWB.SpikeTimes{:}); spk_events = arrayfun(@(r) ...
    cat(2, spk_events{r,:}), (1:size(spk_events,1))', fnOpts{:});
% Mean waveforms
% mu_waves = cellfun(@(w) mean(w, 2, "omitnan"), clWaveforms(:,2), fnOpts{:});
% mu_waves = cat(2, mu_waves{:});
% Main channel
cl_colNames = clInfo.Properties.VariableNames; 
chFlag = contains(cl_colNames, ["channel", "ch"]);
channels = clInfo{:, chFlag};
% Session references
ses_name = arrayfun(@(id) sprintf("unit%d", id), ids);
ses_ref = arrayfun(@(id) types.untyped.ObjectView(['/analysis/', char(id)]), ...
    ses_name);
%% Assigning unit info to NWB
nwbObj.units = types.core.Units('colnames',...
    {'spike_times', 'trials'},...
    'description', 'Trigger-relative Spike Events');
nwbObj.units.spike_times = types.hdmf_common.VectorData(...
    'description', 'Trigger-relative timestamps of spikes');

for cu=1:length(ids)
    eventTrials = tIDs{cu};
    eventTimes = spk_events{cu};
    % waveform = mu_waves(:,cu);
    channel = channels(cu);
    
    % add waveform data to "unitx" and associate with "waveform" column as ObjectView.
    ses = types.core.SpikeEventSeries(...
        'starting_time_rate', fs, ...
        'control', ids(cu),...
        'control_description', 'Units Table ID',...
        'data', eventTimes , ...
        'description', sprintf('Single unit %d', ids(cu)), ...
        'timestamps', eventTimes, ...
        'timestamps_unit', 'seconds',...
        'electrodes', types.hdmf_common.DynamicTableRegion(...
            'description', 'Electrodes involved with these spike events',...
            'table', types.untyped.ObjectView('/general/extracellular_ephys/electrodes'),...
            'data', channel));
    
    nwbObj.analysis.set(char(ses_name(cu)), ses);
    nwbObj.units.addRow('id', ids(cu), 'trials', eventTrials, ...
        'spike_times', eventTimes);
    
    %add this timeseries into the trials table as well. (From tutorial)
    [s_trials, ~, trials_to_data] = unique(eventTrials);
    for j=1:length(s_trials)
        trial = s_trials(j);
        j_loc = j == trials_to_data;
        t_start = find(j_loc, 1);
        t_count = sum(j_loc);
        
        trial_timeseries{trial}(end+1, :) = ...
            {int32(t_start) int32(t_count) ses_ref(cu)};
    end
end
end