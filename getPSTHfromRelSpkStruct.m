function [PSTH, psthTx, Na] = getPSTHfromRelSpkStruct(relSpkStrct, confgStruct)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
catClSpks = @(x, c) cat(2, x.SpikeTimes{c,:});
htOpts = {'BinWidth', confgStruct.BinSize_s,...
    'BinLimits', confgStruct.Viewing_window_s};
fnOpts = {'UniformOutput', false};
Ncond = length(relSpkStrct);
[Ncl, Na] = arrayfun(@(x) size(x.SpikeTimes), relSpkStrct);
timeLapse = confgStruct.Viewing_window_s;
Npsth = round(diff(timeLapse)/confgStruct.BinSize_s);
mdl = fit_poly([1,Npsth], ...
    timeLapse + [1,-1]*(confgStruct.BinSize_s / 2), 1);
psthTx = ((1:Npsth).^[1;0])'*mdl;

PSTH = arrayfun(@(x) arrayfun(@(y) histcounts(catClSpks(relSpkStrct(x), y),...
    htOpts{:}), 1:Ncl(x), fnOpts{:}), 1:Ncond, fnOpts{:});
PSTH = cellfun(@(x) cat(1, x{:}), PSTH, fnOpts{:});
PSTH = cat(3, PSTH{:});
end