function [psth_u, psth_tx, Na] = getPopPSTH_perTrial(rst, cf)
psth_u = [];
if ~checkRelSpkTmsStruct(rst)
    fprintf(1, 'No good relative spike time structure given!\n')
    return 
end
if ~checkConfigStruct(cf)
    fprintf(1, 'No good configuration structure given!\n')
    return
end

binSize = cf.BinSize_s; vw = cf.Viewing_window_s;
fnOpts = {'UniformOutput', false};
histOpts = {'BinWidth', binSize, 'BinLimits', vw};

Ncl = size(rst(1).SpikeTimes, 1);
psth_u = arrayfun(@(c) arrayfun(@(tr) ...
    histcounts([rst(c).SpikeTimes{:,tr}], ...
    histOpts{:}), 1:size(rst(c).SpikeTimes, 2), ...
    fnOpts{:}), 1:length(rst), fnOpts{:});
psth_u = cellfun(@(c) cat(1, c{:})/(Ncl*binSize), psth_u, fnOpts{:});

[Na, Nt] = size(psth_u{1}); Nc = length(psth_u);
mdl_tm_ax = fit_poly([1, Nt], vw + [1, -1]*binSize, 1);
psth_tx = ((1:Nt)'.^[1,0]) * mdl_tm_ax;

end