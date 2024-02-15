function [PSTH_unit_trial, psth_tx, Na] = ...
    getPSTH_perU_perT(relSpkTmsStruct, confStruct, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Emilio Isaías-Camacho @GrohLab 2022
p = inputParser();

addRequired(p, 'relSpkTmsStruct', @checkRelSpkTmsStruct)
addRequired(p, 'confStruct', @checkConfigStruct)
addParameter(p, 'Filter', false, @(x) all([islogical(x), numel(x) == 1]))

parse(p, relSpkTmsStruct, confStruct, varargin{:});

relSpkTmsStruct = p.Results.relSpkTmsStruct;
confStruct = p.Results.confStruct;
fFilt = p.Results.Filter;

timeLapse = confStruct.Viewing_window_s;
binSz = confStruct.BinSize_s;
fnOpts = {'UniformOutput', false};
psthOpts = {'BinLimits', timeLapse, 'BinWidth', binSz, 'Normalization', ...
    'count'};

if fFilt
    % Convolution kernel: Gaussian window
    conv_kernel = gausswin(5, pi);
    conv_kernel = conv_kernel/sum(conv_kernel);

    % Creating a smooth histogram per trial per unit.
    PSTH_unit_trial = arrayfun(@(x) ...
        cellfun(@(y) conv(histcounts(y, psthOpts{:}), conv_kernel, 'same'), ...
        x.SpikeTimes, fnOpts{:}), ...
        relSpkTmsStruct, fnOpts{:});
else
    PSTH_unit_trial = arrayfun(@(x) ...
        cellfun(@(y) histcounts(y, psthOpts{:}), x.SpikeTimes, fnOpts{:}), ...
        relSpkTmsStruct, fnOpts{:});
end

% Converting trials of a single unit in a matrix
PSTH_unit_trial = cellfun(@(x) ...
    arrayfun(@(y) cat(1,x{y,:}), (1:size(x, 1))', fnOpts{:}), ...
    PSTH_unit_trial, fnOpts{:});

% Converting all single unit trial matrix into a 3D matrix:
% Trials x Bins x Units
PSTH_unit_trial = cellfun(@(x) cat(3, x{:}), ...
    PSTH_unit_trial, fnOpts{:});
end