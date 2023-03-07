function [outputArg1,outputArg2] = assignSpkTms2NWB(nwbObj, rstNWB)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Units
nwb.units = types.core.Units('colnames',...
    {'spike_times', 'trials', 'waveforms'},...
    'description', 'Analysed Spike Events');
esHash = data.eventSeriesHash;
ids = regexp(esHash.keyNames, '^unit(\d+)$', 'once', 'tokens');
ids = str2double([ids{:}]);
nwb.units.spike_times = types.hdmf_common.VectorData(...
    'description', 'timestamps of spikes');

for i=1:length(ids)
    esData = esHash.value{i};
    % add trials ID reference
    
    good_trials_mask = ismember(esData.eventTrials, nwb.intervals_trials.id.data);
    eventTrials = esData.eventTrials(good_trials_mask);
    eventTimes = esData.eventTimes(good_trials_mask);
    waveforms = esData.waveforms(good_trials_mask,:);
    channel = esData.channel(good_trials_mask);
    
    % add waveform data to "unitx" and associate with "waveform" column as ObjectView.
    ses = types.core.SpikeEventSeries(...
        'control', ids(i),...
        'control_description', 'Units Table ID',...
        'data', waveforms .', ...
        'description', esHash.descr{i}, ...
        'timestamps', eventTimes, ...
        'timestamps_unit', data.timeUnitNames{data.timeUnitIds(esData.timeUnit)},...
        'electrodes', types.hdmf_common.DynamicTableRegion(...
            'description', 'Electrodes involved with these spike events',...
            'table', types.untyped.ObjectView('/general/extracellular_ephys/electrodes'),...
            'data', channel - 1));
    ses_name = esHash.keyNames{i};
    ses_ref = types.untyped.ObjectView(['/analysis/', ses_name]);
    if ~isempty(esData.cellType)
        ses.comments = ['cellType: ' esData.cellType{1}];
    end
    nwb.analysis.set(ses_name, ses);
    nwb.units.addRow(...
        'id', ids(i), 'trials', eventTrials,'spike_times', eventTimes, 'waveforms', ses_ref);
        %'tablepath', '/units');
    
    %add this timeseries into the trials table as well.
    [s_trials, ~, trials_to_data] = unique(eventTrials);
    for j=1:length(s_trials)
        trial = s_trials(j);
        j_loc = j == trials_to_data;
        t_start = find(j_loc, 1);
        t_count = sum(j_loc);
        
        trial_timeseries{trial}(end+1, :) = {ses_ref int64(t_start) int64(t_count)};
    end
end

end