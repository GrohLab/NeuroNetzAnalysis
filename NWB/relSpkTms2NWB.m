function [outputArg1,outputArg2] = relSpkTms2NWB(relSpkTms, confStruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
% You can find more information about hashes and how they're used on the
% <https://crcns.org/data-sets/motor-cortex/alm-3/about-alm-3 ALM-3 about page>.
fprintf('Processing Data Structure `%s`\n', datastructure_loc);
loaded = load(datastructure_loc, 'obj');
data = loaded.obj;

% wherein each cell is one trial. We must populate this way because trials
% may not be in trial order.
% Trial timeseries will be a compound type under intervals/trials.

%% Units
units = types.core.Units('colnames',...
    {'spike_times', 'trials', 'waveforms'},...
    'description', 'Analysed Spike Events');

units.spike_times = types.hdmf_common.VectorData(...
    'description', 'timestamps of spikes');

% Line in a loop
units.addRow(...
        'id', ids(i), 'trials', eventTrials,'spike_times', eventTimes, 'waveforms', ses_ref,...
        'tablepath', '/units');

trials = types.core.TimeIntervals( ...
    'colnames', {'start_time', 'stop_time'}, ...
    'description', 'trial data and properties', ...
    'id', types.hdmf_common.ElementIdentifiers('data', 0:99), ...
    'start_time', types.hdmf_common.VectorData( ...
        'data', repmat(-0.3,1,100), ...
        'description','start time of trial in seconds' ...
    ), ...
    'stop_time', types.hdmf_common.VectorData( ...
        'data', repmat(0.4,1,100), ...
        'description','end of each trial in seconds' ...
    ));

for ct = 1:Nt

end
end