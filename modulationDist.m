function [frHist, frCents, structString, MI] =...
    modulationDist(fr, frNbin, figTtl, figDir, chExp, structString)
%MODULATIONDIST provides the modulation distribution for the given Nx2
%firing frequency matrix. It returns the histogram for the given bins and
%the center for each bin.
%       [frHist, frCents] = modulationDist(fr, frNbin, figTtl, figDir, *chExp)
%   INTPUS:
%       fr - NxC matrix containing the firing frequency of N clusters with
%            the first column being the control condition and the second
%            being the treatment conditions.
%       frNbin - number of bins for the distribution
%       figTtl - string array specifying the title for the figure and the
%                file name for pdf, emf, and fig files.
%       figDir - path to save the resulting figure.
%       *chExp - 1xM vector representing the chosen experiments
%       *structString
%   OUTPUTS:
%       frHist - 1xfrNbin vector representing the histogram for the
%                modulation index
%       frCents - 1xfrNbin vector representing the bin centers for the
%                 histogram.
%       structString - string as user input
%       MI - NxK array containing the modulation index per combinatorial
%   Modulation Index:
%        $MI = \frac{T - C}{T + C}$
%
% Emilio Isaias-Camacho @GrohLab 2021
%%

fnOpts = {'UniformOutput', 0};

if ~exist('chExp','var')
    chExp = [];
end
% Blue shades for depicting the percentages
qCMap = [224, 231, 250; 123, 152, 234; 21, 50, 132];
qCMap = cat(1, qCMap, flip(qCMap(1:2,:),1))./255;

[Ncl, Ncond] = size(fr);
miCond = nchoosek(1:Ncond,2);

% Modulation index: (Treatment - Control) / (Treatment + Control)
getModulationIndex = @(x) diff(x,1,2)./sum(x,2);
% MI = diff(fr,1,2)./sum(fr,2);
MI = arrayfun(@(x) getModulationIndex(fr(:,miCond(x,:))), 1:size(miCond,1),...
    fnOpts{:}); MI = cat(2, MI{:}); 
sfrFig = figure('Name','Firing rate modulation distribution',...
    'Color',[1,1,1]);
sfrAx = axes('Parent', sfrFig, 'NextPlot', 'add', 'Box', 'off',...
    'Color', 'none');
[frCents, sfrEdg] = prepareLogBinEdges(10.^(-1:1), frNbin);
for cmi = 1:size(MI,2)
    frHist = histcounts(MI(:,cmi), sfrEdg, 'Normalization', 'probability');
    [~, ~, sfrQV, ~] =...  --------------------------------Generalise this part
        exponentialSpread(frHist, frCents, frCents([1,frNbin]));
    plot(sfrAx, frCents, frHist, ':','Color',ones(1,3)*0.75);
    % Depicting the area under 0-5%, 5-25%, 25-75%, 75-95%
    
    qVals10 = [frCents(1), sfrQV([1,2,5,6]), frCents(frNbin)];
    binCents10 = frCents; binCtsQV = interp1(binCents10,...
        frHist, qVals10);
    sdaOpts = {"FaceAlpha", 0.6, "EdgeColor", "none", "FaceColor"};
    for cq = 1:length(qVals10)-1
        qIdx = binCents10 >= qVals10(cq) & binCents10 <= qVals10(cq+1);
        area(sfrAx, [qVals10(cq), binCents10(qIdx), qVals10(cq+1)],...
            [binCtsQV(cq), frHist(qIdx), binCtsQV(cq+1)], sdaOpts{:},...
            qCMap(cq,:));
    end
    % Median , mean, and mode markers
    [~, mxQSub] = max(frHist); mLabels = ["Median","Mean","Mode"];
    triM = [sfrQV(3), mean(MI(:,cmi),'omitnan'), frCents(mxQSub)];
    triBCts = interp1(binCents10, frHist, triM); mCMap = flip(hsv(3),1);
    triMLines = gobjects(size(triM,2),1);
    for cm = 1:length(triM)
        triMLines(cm) = line(sfrAx, repmat(triM(cm),2,1), [0; triBCts(cm)],...
            "Color", mCMap(cm,:), "DisplayName", mLabels(cm)+" "+triM(cm));
    end
    sfrAx.XLim = [-1,1];
    % Axis title and figure file name
    ttlString = figTtl;
    ttlFile = cat(2, ttlString, sprintf(' %d', chExp));
end
if ~exist('structString','var')
    structString = [];
    structAns = inputdlg('What structure are you looking at?', 'Structure');
    if ~isempty(structAns)
        structString = structAns{:};
    end
end
ttlString = cat(2,ttlString, sprintf(' (%s)', structString));
ttlFile = cat(2, ttlFile, sprintf(' (%s)', structString));
lgnd = legend(sfrAx, triMLines(1:length(triM))); 
lgnd.Box = 'off'; lgnd.Location = 'best';
title(sfrAx, ttlString)
xlabel(sfrAx, 'Modulation index (MI)');ylabel(sfrAx, "Population proportion");
saveFigure(sfrFig, fullfile(figDir, ttlFile));
end