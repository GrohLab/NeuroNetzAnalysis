function [figs] = plotLogPSTH(PSTHstruct, varargin)
%PLOTLOGPSTH takes the logPSTH structure generated by getLogTimePSTH and
%plots the PSTH per condition, and their mean.
%   Detailed explanation goes here
%Emilio Isaias-Camacho @ GrohLab 2021

%% Parse inputs
p = inputParser;

% Required arguments auxiliary function and variables
pfn = ["LogPSTH", "Log10TimeAxis", "TimeAxis", "ConditionNames",...
    "DeltaLogStep", "Normalization", "Log10BinEdges"];
checkStruct = @(x) all([isstruct(x), isfield(x, pfn)]);

% Required parameters
% none so far

% Optional parameters
% none so far, maybe ordering or selecting different dimensions from the
% matrix (specific clusters and or conditions?)

p.addRequired('PSTHstruct', checkStruct);

p.KeepUnmatched = true;

p.parse(PSTHstruct, varargin{:});

PSTHstruct = p.Results.PSTHstruct;

%% Auxiliary variables
expSubs = @(x) x(1):x(2);
fnOpts = {'UniformOutput',false};
figOpts = {'Visible','on','Color','w'};
getMI = @(x, d) diff(x,1,d)./sum(x,d,'omitnan');

[Ncl, Nbin, Ncond] = size(PSTHstruct.LogPSTH);
tmWinMS = PSTHstruct.TimeAxis([1,Nbin])*1e3;


%% Figure and axes for displaying PSTHs per condition
natFig = figure(figOpts{:}); natAx = gobjects(Ncond+1,1);
natP = {'Parent', natFig};
% Plotting the mean PSTH for all conditions at the bottom of the figure
condPsth = squeeze(mean(PSTHstruct.LogPSTH, 1, 'omitnan')); 
if strcmpi(PSTHstruct.Normalization, 'prob')
    % Probability
    condPsth = condPsth./sum(condPsth, 'omitnan');
%     cbLabel = 'Probability';
else
    % Firing rate with inhomogenous bin sizes
    [~, lgEdg] = prepareLogBinEdges(PSTHstruct.TimeAxis([1,Nbin]), Nbin);
    tmBinWdth = (diff(10.^lgEdg(:)));
    condPsth = condPsth./tmBinWdth;
%     PSTHstruct.LogPSTH = PSTHstruct.LogPSTH.*reshape(condPsth, [1, Nbin, Ncond]);
%     PSTHstruct.LogPSTH = PSTHstruct.LogPSTH./tmBinWdth';
%     cbLabel = 'Firing rate [Hz]';
end
logMeanEdges = [2,1; 3,0]*[Ncond;1];
natAx(Ncond + 1) = subplot(3, Ncond, expSubs(logMeanEdges), natP{:});
semilogx(natAx(Ncond + 1), PSTHstruct.TimeAxis*1e3, condPsth);
xticklabels(natAx(Ncond + 1), xticks(natAx(Ncond + 1)));
xlim(natAx(Ncond + 1), tmWinMS); xlabel(natAx(Ncond + 1), 'Log time [ms]');
if strcmpi(PSTHstruct.Normalization, 'prob')
    ylabel(natAx(Ncond + 1), 'Firing probability [p(spike|bin)]')
    yyaxis(natAx(Ncond + 1), 'right');
    natAx(Ncond + 1).ColorOrder = lines(Ncond);
    natAx(Ncond + 1).LineStyleOrder = '--';
    semilogx(natAx(Ncond + 1), PSTHstruct.TimeAxis*1e3, cumsum(condPsth));
    natAx(Ncond + 1).YAxis(2).Color = [0.3,0.3,0.3];
else
    ylabel(natAx(Ncond + 1), 'Firing rate [Hz]')
end
xtks = xticks(natAx(Ncond + 1)); 
lgnd = legend(natAx(Ncond + 1), PSTHstruct.ConditionNames.cellstr);
set(lgnd,'Location','best','Box','off');
box(natAx(Ncond + 1), 'off'); 
% Plotting PSTH per cluster
tx = PSTHstruct.Log10TimeAxis; mxClr = max(PSTHstruct.LogPSTH(:));
clrMp = rocket(2^8);
for ccond = 1:Ncond
    imgMat = [ccond, 0; ccond, Ncond];
    natAx(ccond) = subplot(3, Ncond, imgMat * [1;1], natP{:}); 
    imagesc(natAx(ccond), 'XData', tx, 'YData', 1:Ncl,...
        'CData', PSTHstruct.LogPSTH(:,:,ccond),[0, mxClr]);
    colormap(natAx(ccond), clrMp)
    xticklabels(natAx(ccond), 10.^(xticks(natAx(ccond))+3));
    title(natAx(ccond), PSTHstruct.ConditionNames(ccond));
    ylim(natAx(ccond), [1, Ncl]); box(natAx(ccond), 'off');
    xticks(natAx(ccond), log10(xtks)-3); 
    xticklabels(natAx(ccond), 10.^(xticks+3))
end
ylabel(natAx(1), 'Clusters'); 
% cb = colorbar(natAx(ccond), 'Location', 'west', 'Color', [0.8,0.8,0.8],...
%     'Limits', [0, max(condPsth(:))]);
% cb.Label.String = cbLabel;
arrayfun(@(x) xlim(x, PSTHstruct.Log10TimeAxis([1,Nbin])), natAx(1:Ncond));
arrayfun(@(x) set(x.YAxis, 'Visible', 'off'),...
    natAx(setdiff(1:(Ncond + 1), [1, Ncond + 1])), fnOpts{:});
arrayfun(@(x) set(x,'Color','none'), natAx);
natFig.Visible = 'on';
figs(1) = natFig;
%% Figure for displaying a comparison between condition permutations
if Ncond > 1
    posBar = {'FaceColor', 'g', 'EdgeColor', 'none'};
    negBar = {'FaceColor', 'r', 'EdgeColor', 'none'};
    permFig = figure(figOpts{:}); permP = {'Parent', permFig};
    % Comparison between conditions
    permCond = nchoosek(1:Ncond,2); Nperm = size(permCond,1);
    prmAx = gobjects(Nperm*2,1);
    for cperm = 1:Nperm
        % MI per cluster
        imgMat = [cperm, 0; cperm, Nperm];
        prmAx(cperm) = subplot(3, Nperm, imgMat * [1;1], permP{:});
        cMI = getMI(PSTHstruct.LogPSTH(:,:,permCond(cperm,:)),3);
        cMI(isnan(cMI)) = 0; 
        imagesc(prmAx(cperm), 'XData', tx, 'YData', 1:Ncl, 'CData', cMI);
        colormap(prmAx(cperm), traffic(101)); ylim(prmAx(cperm), [1,Ncl]);
        title(prmAx(cperm), sprintf('%s vs %s', ...
            PSTHstruct.ConditionNames(permCond(cperm,[2,1]))))
        % MI per condition all clusters
        miSpSub = Nperm + cperm;
        prmAx(miSpSub) = subplot(3, Nperm, [2, cperm]*[Nperm;1], permP{:});
        condMI = getMI(condPsth(:,permCond(cperm,:)),2);
        condMI(isnan(condMI)) = 0;
        bar(prmAx(miSpSub), tx(condMI>0), condMI(condMI>0), posBar{:});
        set(prmAx(miSpSub),'NextPlot','add');
        bar(prmAx(miSpSub), tx(condMI<=0), condMI(condMI<=0), negBar{:});
        set(prmAx(miSpSub),'Color','none');
        xticks(prmAx(miSpSub), log10(xtks)-3); 
        xticklabels(prmAx(miSpSub), 10.^(xticks+3))
    end
    ylabel(prmAx(1), 'Clusters')
    arrayfun(@(x) xlim(x, PSTHstruct.Log10TimeAxis([1,Nbin])), prmAx(1:Nperm*2));
    arrayfun(@(x) set(x.XAxis,'Visible','off'), prmAx(1:Nperm));
    arrayfun(@(x) set(x.YAxis,'Visible','off'), prmAx(2:Nperm));
    arrayfun(@(x) box(x, 'off'), prmAx);
    linkaxes(prmAx, 'x'); linkaxes(prmAx(Nperm+1:Nperm*2), 'xy')
    lgnd = legend(prmAx(Nperm*2), {'Potentiation','Depression'});
    lgnd.Location = 'best'; lgnd.Box = 'off';
    ylabel(prmAx(Nperm+1), 'Modulation Index')
    permFig.Visible = 'on';
    figs(2) = permFig;
end
end

