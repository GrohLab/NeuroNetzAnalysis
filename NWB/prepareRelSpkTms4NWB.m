function [rstNWB] = prepareRelSpkTms4NWB(sessionPath, configStruct, varargin)
%prepareRelSpkTms4NWB prepares the spike times in the structure to add them
%to the NWB framework.
%   [rstNWB] = prepareRelSpkTms4NWB(sessionPath, configStruct, varargin)

%% Input validation
p = inputParser();
necessaryFields = {'Experiment', 'Viewing_window_s', 'Response_window_s', ...
    'Spontaneous_window_s', 'BinSize_s', 'Trigger', 'ConsideredConditions'};
validateCS = @(x) isstruct(x) & all(isfield(x, necessaryFields));
istxt = @(x) isstring(x) | ischar(x) | iscellstr(x);

addRequired(p, "sessionPath", @(x) exist(x, "dir"))
addRequired(p, "configStruct", validateCS)
addParameter(p, 'ConditionSelect', "all", istxt)
parse(p, sessionPath, configStruct, varargin{:});

sessionPath = p.Results.sessionPath;
configStruct = p.Results.configStruct;
condSel = p.Results.ConditionSelect;

%% Auxiliary functions and variables
fnOpts = {'UniformOutput', false};
pathHere = @(x) fullfile(sessionPath, x);
look4this = @(x) dir(pathHere(x));
expandName = @(x) fullfile(x.folder, x.name);

[~, baseName] = fileparts(configStruct.Experiment);
rstPttrn = string(baseName)+ " " + sprintf("RW%.2f - %.2f ms ", ...
    configStruct.Response_window_s.*1e3) + sprintf("SW%.2f - %.2f ms ", ...
    configStruct.Spontaneous_window_s.*1e3) + sprintf("VW%.2f - %.2f ms ", ...
    configStruct.Viewing_window_s.*1e3) + string(configStruct.Trigger.Name) + ...
    " (unfiltered) exportSpkTms.mat";

rstFiles = look4this(rstPttrn); ciFiles = look4this("cluster_info.tsv");
if any(~arrayfun(@(x) exist(expandName(x), "file"), [rstFiles; ciFiles]))
    fprintf(1, "Unable to proceed with these configuration parameters"+...
        " or session path!\n")
    return
end
%% Loading and processing data
clInfo = getClusterInfo(expandName(ciFiles));
info_col_names = clInfo.Properties.VariableNames;
idFlag = ismember(info_col_names, ["cluster_id", "id"]);
auFlag = ismember(info_col_names, "ActiveUnit"); 
gclID = clInfo{clInfo.(info_col_names{auFlag})==1,info_col_names{idFlag}};

load(expandName(rstFiles), "relativeSpkTmsStruct")
condNames = arrayfun(@(x) string(x.name), relativeSpkTmsStruct);
if strcmp(condSel, "all")
    % All conditions into the NWB file
    condSel = condNames;
else
    % Selected conditions by their name or it might be a cell array
    % with char arrays or strings
    condSel = condSel(ismember(condNames, condSel));
end

condFlags = arrayfun(@(rst) ismember(string(rst.name), string(condSel)), ...
    relativeSpkTmsStruct);
[Ncl, Na] = arrayfun(@(rst) size(rst.SpikeTimes), relativeSpkTmsStruct(condFlags));
spk_per_trial_per_unit = arrayfun(@(rst) cellfun(@(t) numel(t), ...
    rst.SpikeTimes), relativeSpkTmsStruct(condFlags), fnOpts{:});
offst = [0, cumsum(Na(1:end-1))]; 
% Assigning a trial ID per spike in the structure.
trialID_per_unit = cellfun(@(stu, si) arrayfun(@(r) arrayfun(@(rp,v) ...
    repmat(v,1,rp), stu(r,:), (1:size(stu,2)) + offst(si), fnOpts{:}), ...
    1:size(stu,1), fnOpts{:}), spk_per_trial_per_unit, ...
    num2cell(1:size(spk_per_trial_per_unit,2)), fnOpts{:});
trialID_per_unit = cellfun(@(x) cat(1, x{:}), trialID_per_unit, fnOpts{:});
trialID_per_unit = cellfun(@(stu) arrayfun(@(x) cat(2, stu{x,:}), ...
    (1:size(stu,1))', fnOpts{:}), trialID_per_unit, fnOpts{:});
% Linearising the spike times for all units and all trials, per condition
lin_spkTms = arrayfun(@(rst) arrayfun(@(t) cat(2, rst.SpikeTimes{t,:}), ...
    (1:size(rst.SpikeTimes,1))', fnOpts{:}), relativeSpkTmsStruct(condFlags), ...
    fnOpts{:});
% NWB output
rstNWB = struct('TrialType', string(condSel), 'TrialID', {trialID_per_unit},...
    'UnitID',gclID,'SpikeTimes', {lin_spkTms}, 'NumTrials', Na, 'NumUnits', Ncl(1));
end