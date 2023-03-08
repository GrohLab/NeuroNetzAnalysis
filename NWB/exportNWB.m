function [nwbObj] = exportNWB(nwbObj, trial_timeseries, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Input parsing
p = inputParser;

istxt = @(x) isstring(x) | ischar(x);

addRequired(p, "nwbObj", @(x) isa(x, "NwbFile"))
addRequired(p, "trial_timeseries", @iscell)
addParameter(p, "OutDir", "", istxt)

parse(p, nwbObj, trial_timeseries, varargin{:})

nwbObj = p.Results.nwbObj;
trial_timeseries = p.Results.trial_timeseries;
outDir = p.Results.OutDir;
%% 
fnOpts = {'UniformOutput', false};
if ~exist(outDir,"dir")
    answ = questdlg("Path doesn't exist. Would you like to create such folder?",...
        "Out directory", "Yes", "No", "Yes");
    if strcmpi(answ,"yes")
        if ~mkdir(outDir)
            fprintf(1, "Couldn't create folder!\n")
            outDir = pwd;
        end
    else
        fprintf(1, "Exiting\n")
        return
    end
end

%first, we'll format and store |trial_timeseries| into |intervals_trials|.
% note that |timeseries_index| data is 0-indexed.
ts_len = cellfun("size", trial_timeseries, 1);
is_len_nonzero = ts_len > 0;
ts_len_nonzero = ts_len(is_len_nonzero);
nwbObj.intervals_trials.timeseries_index.data = cumsum(ts_len_nonzero);
% intervals/trials/timeseries is a compound type so we use cell2table to
% convert this 2-d cell array into a compatible table.
trial_timeseries = cat(1, trial_timeseries{is_len_nonzero});
trial_timeseries = arrayfun(@(c) cat(1, trial_timeseries{:,c}), ...
    1:size(trial_timeseries, 2), fnOpts{:});

nwbObj.intervals_trials.timeseries.data = table(trial_timeseries{:}, ...
    'VariableNames', {'idx_start', 'count', 'timeseries'});

outDest = fullfile(outDir, [nwbObj.identifier '.nwb']);
if 2 == exist(outDest, 'file')
    delete(outDest);
end
nwbExport(nwbObj, cellstr(outDest));
end