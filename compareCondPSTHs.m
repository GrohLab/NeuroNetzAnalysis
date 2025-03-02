function [ppFig, PSTHall] = compareCondPSTHs(PSTH, Na, binSz, vWin, condNames)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% fnOpts = {'UniformOutput', false};
axOpts = {'Box', 'off', 'Color', 'none'};
lgOpts = [axOpts(:)', {'Location'}, {'best'}, {'AutoUpdate'}, {'off'}];
zlOpts = {'LineStyle','--','Color'};
rwN = 4; axSbs = (0:rwN-2)'; 
[Nt, Nccond] = size(PSTH, [2,3]);
psthTx = (0:Nt-1) * binSz + vWin(1);
PSTHall = PSTH ./ reshape(Na*binSz, 1, 1, Nccond);
zPSTH = zscore(PSTHall, 1, [2,3]); zPopPSTH = squeeze(mean(zPSTH, 1));
if isvector(zPopPSTH) && isrow(zPopPSTH)
    zPopPSTH = zPopPSTH(:);
end
Mxe = max(zPopPSTH, [], "all"); Mne = min(zPopPSTH, [], "all");
figOpts = {'Visible', 'on'};
if ~strcmp( computer, 'PCWIN64' )
    figOpts(2) = {'off'};
end
ppFig = figure('Color', 'w', 'Name', 'Conditions PSTH', figOpts{:} );
axs = gobjects(Nccond+1,1); 

try
    % Decreasing the contrast of colormap since the highest value is a
    % bright yellow that is almost invisible
    lnsCM = plasma( ceil( Nccond *1.1));
    clrMap = rocket;
catch
    fprintf(1, "tsipkens/cmap repository not in MATLAB's path!\n")
    fprintf(1, "Using copper colormap.\n")
    lnsCM = lines(Nccond);
    clrMap = hot;
end


for cc = 1:Nccond
    axs(cc) = subplot(rwN, Nccond, (axSbs.^[1,0])*[Nccond;cc], ...
        'Parent', ppFig);
    imagesc(axs(cc), vWin*1e3, [], zPSTH(:,:,cc), [Mne, Mxe]);
    xlabel(axs(cc), 'Time [ms]'); yticks(axs(cc), []);
    title(axs(cc), "{\color[rgb]{"+join(string(lnsCM(cc,:)))+"}"+ ...
        condNames{cc}+"}"); colormap(clrMap);
    xline(0, zlOpts{:}, 0.75*ones(1,3))
end
ylabel(axs(1), 'Units');
axs(Nccond+1) = subplot(rwN, Nccond, 1+Nccond*(rwN - 1):Nccond*rwN, ...
    'NextPlot', 'add');
lObj = arrayfun(@(x) plot(axs(end), 1e3*psthTx, zPopPSTH(:,x), ...
    "Color", lnsCM(x,:), "LineWidth", 1.5, ...
    "DisplayName", condNames{x}), 1:Nccond);
lgnd = legend(axs(end), lObj); set(lgnd, lgOpts{:})
xline(0, zlOpts{:}, 0.35*ones(1,3))
xlabel(axs(end), 'Time [ms]'); ylabel(axs(end), "Z-score")
arrayfun(@(x) set(get(x, 'YAxis'), 'Visible','off'), ...
    axs(setdiff(1:(Nccond + 1), [1, Nccond + 1])));
set(axs, axOpts{:}); linkaxes(axs, 'x'); xlim(axs, vWin*1e3)
axPos = axs(Nccond).Position; cbPos = [0.9589, axPos(2), 0.0200, axPos(4)];
cb = colorbar(axs(Nccond), axOpts{1:2}, 'Position', cbPos, ...
    'TickDirection', 'none', 'AxisLocation', 'in', 'Color', 0.15*ones(1,3));
cb.Label.String = 'Z-score';
end