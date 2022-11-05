function [fig, FRpu] = plotISI(spkSubs, fs, ID)
% PLOTISI takes the spikes subscripts and displays the inter-spike
% intervals in a continuous histogram y-log(x) manner.
%   Input parameters:
%       - spkSubs - cell array containing the spike subscripts for each
%       unit
%       - fs - Sampling frequency used to acquire the signals
%   Output parameters:
%       - fig - figures array containing all created figures.
%       - FRpu - Firing rate per unit
%% Input validation
fig = gobjects(1);
[Nr, Nc] = size(spkSubs);
if iscell(spkSubs) && all(cellfun(@isnumeric,spkSubs))
    Nu = max(Nr,Nc);
    % Validation for subscripts or time points.
    SubsFlag = any(cellfun(@all,cellfun(@eq,...
        cellfun(@minus,...
        cellfun(@round,spkSubs,...
        'UniformOutput',false),spkSubs,...
        'UniformOutput',false),repmat({0},Nr,Nc),...
        'UniformOutput',false)));
    mintd = 1/fs;
    if ~SubsFlag
        fprintf(1,'Given time points instead of subscripts.')
        fprintf(1,' Considering time in seconds\n')
        fs = 1;
    end
else
    fprintf(1,'Not sure how to manage this entry...\n')
    fprintf(1,'Please provide either a cell array, a numeric vector, or a ')
    fprintf(1,'logical time series\n')
    return
end
IDflag = false;
if exist('ID','var') && Nu == numel(ID)
    fprintf(1,'Naming each unit appropiately\n')
    IDflag = true;
end

fig = gobjects(Nu,1);
if Nu >= 100
    warning('More than or equal to 100 units. Consider curation')
end
%% Inter-spike intervals
FRpu = zeros(Nu,3);
for cu = 1:Nu
    lisi = log10(diff(spkSubs{cu}./fs));
    if any(isinf(lisi)) 
        fprintf(1,'Cluster %d ',cu);
        if IDflag
            fprintf(1,'(%s)',ID{cu})
        end
        fprintf(1,' has repeated time points.\n')
        fprintf(1,'These time points will not be considered.\n')
        lisi(isinf(lisi)) = [];
    end
    if isempty(lisi)
        fprintf(1,'Cluster %d ',cu);
        if IDflag
            fprintf(1,'(%s)',ID{cu})
        end
        fprintf(1,' has one time point.\n')
        fprintf(1,'Skipping...\n')
        continue
    end
    FRpu(cu,:) = [1/(10^mean(lisi,'omitnan')),...
        1/(10^median(lisi,'omitnan')),...
        1/(10^mode(lisi))];
    figure('Visible','off','Color',[1,1,1]);
    hisi = histogram(lisi,linspace(log10(mintd),3,100));
    cts = hisi.BinCounts;
    bns = (hisi.BinEdges(1:end-1) + hisi.BinEdges(2:end))/2;
    plot(bns,cts./sum(cts),'LineWidth',1);
    fig(cu) = gcf;
    ax = fig(cu).Children;
    ax.XTickLabel = 10.^cellfun(@str2double,ax.XTickLabel) * 1e3;
    xlabel(ax,'Time [ms]'); ylabel(ax,'ISI Probability');
    if IDflag
        title(ax,['Unit #',num2str(cu),' (',ID{cu},')']);
    else
        title(ax,['Unit #',num2str(cu)]);
    end
    grid(ax,'on')
    Ncts = cts/sum(cts);
    yyaxis('right');plot(bns,cumsum(Ncts),'LineStyle','--')
    ylabel('Cumulative fraction');ax = fig(cu).Children;
    ax.YAxis(2).Limits = [0, 1];
    ax.YAxis(2).Color = [0.1,0.1,0.1];
    fig(cu).Visible = 'on';
end
end