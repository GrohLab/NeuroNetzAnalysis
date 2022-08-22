function [PSTH_unit_trial] = getPSTH_perU_perT(relSpkTmsStruct, confStruct)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

timeLapse = confStruct.Viewing_window_s;
binSz = confStruct.BinSize_s;

fnOpts = {'UniformOutput', false};
psthOpts = {'BinLimits', timeLapse, 'BinWidth', binSz, 'Normalization', ...
    'countdensity'};

% Convolution kernel: Gaussian window
conv_kernel = gausswin(5, pi);
conv_kernel = conv_kernel/sum(conv_kernel);

% Creating a smooth histogram per trial per unit.
PSTH_unit_trial = arrayfun(@(x) ...
    cellfun(@(y) conv(histcounts(y, psthOpts{:}), conv_kernel, 'same'), ...
    x.SpikeTimes, fnOpts{:}), ...
    relSpkTmsStruct, fnOpts{:});

% Converting trials of a single unit in a matrix
PSTH_unit_trial = cellfun(@(x) ...
    arrayfun(@(y) cat(1,x{y,:}), (1:size(x, 1))', fnOpts{:}), ...
    PSTH_unit_trial, fnOpts{:});

% Converting all single unit trial matrix into a 3D matrix:
% Trials x Bins x Units
PSTH_unit_trial = cellfun(@(x) cat(3, x{:}), ...
    PSTH_unit_trial, fnOpts{:});
end