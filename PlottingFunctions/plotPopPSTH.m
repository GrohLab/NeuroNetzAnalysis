function fig_out = plotPopPSTH(PSTHStruct, varargin)
%PLOTPOPPSTH This function accepts the PSTH structure created by getPopPSTH
%and plots all PSTHs of all experiments together with a population PSTH at
%the bottom of the figure. It can plot the PSTH structure in a given axis
%or create a new figure to plot this results.
%   Emilio Isa??as-Camacho @ GrohLab 2019

%% Input validation
if isstruct(PSTHStruct) && isfield(PSTHStruct,'TimeAxis')
    p = inputParser;
    % Bin size validation
    defaultBin = PSTHStruct.BinSize;
    small_limit = mean(diff(PSTHStruct.TimeAxis));
    attributes = {'scalar','>',small_limit};%,...
        ...'<=', sum(abs([PSTHStruct.TimeAxis(1),PSTHStruct.TimeAxis(end)]))*0.1};
    addOptional(p,'BinSize', defaultBin, @(x)validateattributes(x,...
        {'numeric'},attributes));
    % Structure validation
    addRequired(p,'PSTHStruct',@isstruct);
    % Display (figure)
    figchk = 'matlab.ui.Figure';
    defaultFig = figure('Visible','off',...
        'Color',[1,1,1],...
        'Name','Population PSTH');
    addOptional(p,'Figure',defaultFig,@(x)isa(x,figchk));
    % Plotting style
    defaultStyle = 'simple';
    validStyles = {'simple','middle','crowded'};
    checkStyle = @(x) any(validatestring(x,validStyles));
    addOptional(p,'PlotStyle',defaultStyle,checkStyle);
    % Actual validation
    p.parse(PSTHStruct, varargin{:})
else
    fprintf('The PSTH structure is either corrupt or hasn''t been created\n')
    return
end
%% Initialization
% Bin size
binsize = p.Results.BinSize;
% Figure
fig = p.Results.Figure;
% Plotting style
plotStyle = p.Results.PlotStyle;
%% Plotting population PSTHs
switch plotStyle
    case 'simple'
        fig_out = plotSimple(PSTHStruct, binsize, fig);
    case 'middle'
        fig_out = plotMiddle(PSTHStruct, binsize, fig);
    case 'crowded'
        fig_out = plotCrowded(PSTHStruct, binsize, fig);
end

end

function fig_out = plotSimple(PSTHStruct, binsize, fig)
% Number of (user-defined) conditions combinations
Ncc = length(PSTHStruct.PSTHs);
% Figure organization
st = sqrt(Ncc);
if st - fix(st) == 0
    Nrow = st;
    Ncol = Nrow;
elseif ~mod(Ncc,2) && Ncc/2 <= 5
    Nrow = 2;
    Ncol = Ncc/Nrow;
elseif ~mod(Ncc,3) && Ncc/3 <= 5
    Nrow = 3;
    Ncol = Ncc/Nrow;
else
    Ncol = 5;
    Nrow = ceil(Ncc/Ncol);
end
ttx = PSTHStruct.TimeAxis;
btx = ttx(1):binsize:ttx(end);
[Nv, Ns, Nx] = size(PSTHStruct.PSTHs(1).PSTHstack);
ax = gobjects(Ncc,1);
fs = 1/mean(diff(ttx));
cmap = jet(Nv-2);
% Gray, olive green and blue shades
cmap = [1/3, 0.502, 0;0,0,0;cmap];
fig.Visible = 'on';
for ccc = 1:Ncc
    % Computing the necessary counters and subscripts to manage the
    % multiple subplots in the graphs.
    auxPSTH = zeros(Nv, Ns, 'single');
    ax(ccc) = subplot(Nrow,Ncol,ccc);
    ax(ccc).Box = 'off';
    ax(ccc).Parent = fig;
    ax(ccc).NextPlot = 'add';
    % Condition name
    condNames = PSTHStruct.PSTHs(ccc).ConditionCombination;
    condtitle = condNames{1};
    cv = 2;
    while cv <= numel(condNames)
        condtitle = cat(2,condtitle,' + ',condNames{cv});
        cv = cv + 1;
    end
    tpe = Nx;
    for cexp = 1:Nx
        if PSTHStruct.PSTHs(ccc).TrialsPerExperiment(cexp)
            auxPSTH = auxPSTH + PSTHStruct.PSTHs(ccc).PSTHstack(:,:,cexp)./...
                PSTHStruct.PSTHs(ccc).TrialsPerExperiment(cexp);
        else
            tpe = tpe - 1;
        end
    end
    bPSTH = binTime(auxPSTH(2:Nv,:), binsize, fs);
    % Removing those variables which didn't occure at all
    emptyFlag = sum(bPSTH,1) == 0;
    bPSTH = bPSTH(:,~emptyFlag);
    nNv = size(bPSTH,2);
    lbs = PSTHStruct.SignalIDs([false,~emptyFlag]);
    cSubs = find([false,~emptyFlag]);
    for cvar = 1:nNv
        if any(strcmpi(lbs{cvar},{'spikes','neuron','unit'}))
            yyaxis(ax(ccc),'right');
            bar(ax(ccc),btx,bPSTH(:,cvar)./binsize,...
                'EdgeColor','none','FaceColor',cmap(2,:),...
                'FaceAlpha',0.3,'DisplayName',lbs{cvar},'BarWidth',1);
            ax(ccc).YAxis(2).Limits =...
                [0,max(bPSTH(:,cvar),[],'omitnan')/binsize];
            ax(ccc).YAxis(2).Label.String = 'Hz';
            yyaxis(ax(ccc),'left');
        else
            bar(ax(ccc),btx,bPSTH(:,cvar)./tpe,...
                'EdgeColor','none','FaceAlpha',0.3,...
                'FaceColor',cmap(cSubs(cvar),:),...
                'DisplayName',lbs{cvar});
        end
    end
    plot(ax(ccc),ttx, auxPSTH(1,:)./tpe, 'Color',cmap(1,:),...
                'DisplayName',PSTHStruct.Trigger,'LineStyle','-');
    ax(ccc).Title.String = condtitle;
    
    Nt = sum(PSTHStruct.PSTHs(ccc).TrialsPerExperiment);
    ax(ccc).YAxis(1).Color = [0,0,0];
    try
        ax(ccc).YAxis(2).Color = [0.3,0.3,0.3];
    catch
        fprintf('This experiment is apparently empty.\n')
    end
    ax(ccc).YAxis(1).TickLabels = [];
    ax(ccc).YAxis(1).TickValues = [];
    ax(ccc).YAxis(1).Label.String = sprintf('Experiments: %d/_{trials: %d}',...
        tpe,Nt);
    xlabel(ax(ccc),sprintf('Time_{%.2f ms} [s]',binsize*1e3))
    legend(ax(ccc),'show','Location','best')
    linkaxes(ax,'x')
end
fig_out = fig;
end

function fig_out = plotMiddle(PSTHStruct, binsize, fig)
Ncc = length(PSTHStruct.PSTHs);
if Ncc == 1
    figs = fig;
else
    figs = gobjects(Ncc,1);
    figs(1) = fig;
end
ttx = PSTHStruct.TimeAxis;
btx = ttx(1):binsize:ttx(end);
[Nv, Ns, Nx] = size(PSTHStruct.PSTHs(1).PSTHstack);
Nx_g = Nx + 1;
fs = 1/mean(diff(ttx));
cmap = jet(Nv-2);
% Green, black and jet shades
cmap = [1/3, 0.502, 0;0,0,0;cmap];
for ccc = 1:Ncc
    % Computing the necessary counters and subscripts to manage the
    % multiple subplots in the graphs.
    if ~(figs(ccc).isvalid && isa(figs(ccc),'matlab.ui.Figure'))
        figs(ccc) = figure('Name','Population PSTH',...
            'Color',[1,1,1],'Visible','off');
    end
    auxPSTH = zeros(Nv,Ns,'single');
    ax = gobjects(Nx_g,1);
    Nt = sum(PSTHStruct.PSTHs(ccc).TrialsPerExperiment);
    % Creating the title per condition
    condNames = PSTHStruct.PSTHs(ccc).ConditionCombination;
    condtitle = condNames{1};
    cv = 2;
    while cv <= numel(condNames)
        condtitle = cat(2,condtitle,' + ',condNames{cv});
        cv = cv + 1;
    end
    Nxt = Nx;
    for cexp = 1:Nx_g
        ax(cexp) = subplot(Nx_g,1,cexp,'Parent',figs(ccc),'Box','off');
        % Setting important axis properties to force a 'clean' display
        hold(ax(cexp),'on')
        % Binning the signals using the _binsize_ variable. When the
        % experiments are over, the overal PSTH is computed and displayed
        % using the catch sections.
        
        try
            if logical(PSTHStruct.PSTHs(ccc).TrialsPerExperiment(cexp))
                auxPSTH = auxPSTH +...
                    PSTHStruct.PSTHs(ccc).PSTHstack(:,:,cexp)./...
                    PSTHStruct.PSTHs(ccc).TrialsPerExperiment(cexp);
            else
                Nxt = Nxt - 1;
            end
            [bPSTH, binWidth] =...
                binTime(PSTHStruct.PSTHs(ccc).PSTHstack(2:Nv,:,cexp),...
                binsize,fs);
        catch
            bPSTH = binTime(auxPSTH(2:Nv,:), binsize, fs);
        end
        % Removing those variables which didn't occure at all
        emptyFlag = sum(bPSTH,1) == 0;
        bPSTH = bPSTH(:,~emptyFlag);
        nNv = size(bPSTH,2);
        lbs = PSTHStruct.SignalIDs([false,~emptyFlag]);
        cSubs = find([false,~emptyFlag]);
        % Normalization value (For each experiment and overall in the catch
        % section
        try
            tpe = PSTHStruct.PSTHs(ccc).TrialsPerExperiment(cexp);
        catch
            tpe = Nxt;
        end
        for cvar = 1:nNv
            % bln = repmat(ofst(cexp),1,Nb);
            if any(strcmpi(lbs{cvar},{'spikes','neuron','unit'}))
                yyaxis(ax(cexp),'right');
                bar(ax(cexp),btx,bPSTH(:,cvar)./binsize,...
                    'EdgeColor','none','FaceColor',cmap(2,:),...
                    'FaceAlpha',0.3,'DisplayName',lbs{cvar},'BarWidth',1);
                ax(cexp).YAxis(2).Limits =...
                    [0,max(bPSTH(:,cvar),[],'omitnan')/binsize];
                ax(cexp).YAxis(2).Label.String = 'Hz';
                yyaxis(ax(cexp),'left');
            else                
                bar(ax(cexp),btx,bPSTH(:,cvar)./tpe,...
                    'EdgeColor','none','FaceAlpha',0.3,...
                    'FaceColor',cmap(cSubs(cvar),:),...
                    'DisplayName',lbs{cvar},'BarWidth',1);
            end
        end
        % Display of the unbinned trigger for each experiment and overall
        % in the catch section.
        try
            plot(ax(cexp),ttx,...
                PSTHStruct.PSTHs(ccc).PSTHstack(1,:,cexp)./tpe,...
                'Color',cmap(1,:),'DisplayName',PSTHStruct.Trigger,...
                'LineStyle','-');
        catch
            plot(ax(cexp),ttx, auxPSTH(1,:)./tpe, 'Color',cmap(1,:),...
                'DisplayName',PSTHStruct.Trigger,'LineStyle','-');
        end
        ax(cexp).XAxis.Visible = 'off';
        ax(cexp).YAxis(1).Color = [0,0,0];
        try
            ax(cexp).YAxis(2).Color = [0.3,0.3,0.3];
        catch
            fprintf('This experiment is apparently empty.\n')
        end
        ax(cexp).YAxis(1).TickLabels = [];
        ax(cexp).YAxis(1).TickValues = [];
        ax(cexp).YAxis(1).Label.String = sprintf('%s_{%d}',num2str(cexp),...
            tpe);
        ax(cexp).YAxis(1).Label.Rotation = 0;
    end
    ax(1).Title.String = condtitle;
    ax(Nx_g).YAxis(1).Label.String = sprintf('P_{%d}',Nt);
    ax(Nx_g).XAxis.Visible = 'on';
    legend(ax(Nx_g),'show','Location','bestoutside')
    xlabel(ax(Nx_g),sprintf('Time_{%.2f ms} [s]',binsize))
    linkaxes(ax,'x')
    figs(ccc).Visible = 'on';
end
fig_out = figs;
end
% The plotCrowded function has probably no much use. It contains several
% bugs that might get fixed as I develop the plotPopPSTH function and
% encapsulate several plotting routines into smaller functions.
function fig_out = plotCrowded(PSTHStruct, binsize, fig)
% Number of (user-defined) conditions combinations
Ncc = length(PSTHStruct.PSTHs);
% Figure organization
st = sqrt(Ncc);
if st - fix(st) == 0
    Nrow = st;
    Ncol = Nrow;
elseif ~mod(Ncc,2) && Ncc/2 <= 5
    Nrow = 2;
    Ncol = Ncc/Nrow;
elseif ~mod(Ncc,3) && Ncc/3 <= 5
    Nrow = 3;
    Ncol = Ncc/Nrow;
else
    Ncol = 5;
    Nrow = ceil(Ncc/Ncol);
end
ttx = PSTHStruct.TimeAxis;
btx = ttx(1):binsize:ttx(end);
[Nv, Ns, Nx] = size(PSTHStruct.PSTHs(1).PSTHstack);
ax = gobjects(Ncc,Nx+1);
fs = 1/mean(diff(ttx));
cmap = jet(Nv-2);
% Gray, olive green and blue shades
cmap = [0.1,0.1,0.1;1/3, 0.502, 0;cmap];
fig.Visible = 'on';
% ofst = (0:Nx-1) * 1.2;
Nb = numel(btx);
for ccc = 1:Ncc
    % Computing the necessary counters and subscripts to manage the
    % multiple subplots in the graphs.
    ccol = mod(ccc-1,Ncol)+1;
    auxPSTH = zeros(Nv-1,Ns,'single');
    trigsig = zeros(1,Ns,'single');
    for cexp = 1:Nx+1
        crow = ((ccc - ccol)/Ncol);
        crowFig = (crow*(Nx+1)) + cexp;
        spIdx = (crowFig-1)*Ncol + ccol;
        ax(ccc,cexp) = subplot(Nrow*(Nx+1),Ncol,spIdx);
        % Setting important axis properties to force a 'clean' display
        ax(ccc,cexp).Box = 'off';
        ax(ccc,cexp).Parent = fig;
        hold(ax(ccc,cexp),'on')
        % Binning the signals using the _binsize_ variable. When the
        % experiments are over, the overal PSTH is computed and displayed
        % using the catch sections.
        try
            auxPSTH = auxPSTH + PSTHStruct.PSTHs(ccc).PSTHstack(2:Nv,:,cexp)./...
                PSTHStruct.PSTHs(ccc).TrialsPerExperiment(cexp);
            trigsig = trigsig + PSTHStruct.PSTHs(ccc).PSTHstack(1,:,cexp)./...
                PSTHStruct.PSTHs(ccc).TrialsPerExperiment(cexp);
            [bPSTH, binWidth] =...
                binTime(PSTHStruct.PSTHs(ccc).PSTHstack(2:Nv,:,cexp),...
                binsize,fs);
        catch
            bPSTH = binTime(auxPSTH, binsize, fs);
        end
        % Removing those variables which didn't occure at all
        emptyFlag = sum(bPSTH,1) == 0;
        bPSTH = bPSTH(:,~emptyFlag);
        nNv = size(bPSTH,2);
        lbs = PSTHStruct.SignalIDs([false,~emptyFlag]);
        cSubs = find([false,~emptyFlag]);
        % Normalization value (For each experiment and overall in the catch
        % section
        try
            tpe = PSTHStruct.PSTHs(ccc).TrialsPerExperiment(cexp);
        catch
            tpe = Nx;
        end
        for cvar = 1:nNv
            % bln = repmat(ofst(cexp),1,Nb);
            if any(strcmpi(lbs{cvar},{'spikes','neuron','unit'}))
                yyaxis(ax(ccc,cexp),'right');
                tmx = [btx(Nb), btx(1), btx]';
                y = [0;0;...ofst(cexp)+flip(bPSTH(:,cvar)./tpe)];
                    flip(bPSTH(:,cvar)./tpe)];
                pts = [tmx,y];
                psh = polyshape(pts,'Simplify',false);
                plot(ax(ccc,cexp),psh,...
                    'EdgeColor','none','FaceColor',cmap(2,:),...
                    'DisplayName',lbs{cvar});
                ax(ccc,cexp).YAxis(2).Limits = [0,max(y,[],'omitnan')];
                yyaxis(ax(ccc,cexp),'left');
%                 patch(ax(ccc,cexp),...
%                     [btx(Nb), btx(1), btx],...
%                     [repmat(ofst(cexp),1,2),...
%                     ofst(cexp)+flip(bPSTH(:,cvar)./tpe)'],cmap(cvar,:))
            else
                % area(ax(ccc,cexp),[btx',btx'],[bln',ofst(cexp)+(bPSTH(:,cvar)./tpe)],...
                bar(ax(ccc,cexp),btx,bPSTH(:,cvar)./tpe,...
                    'EdgeColor','none','FaceAlpha',0.3,...
                    'FaceColor',cmap(cSubs(cvar),:),...
                    'DisplayName',lbs{cvar});%,'BaseValue',ofst(cexp),...
                %'ShowBaseLine','off')
            end
        end
        % Display of the unbinned trigger for each experiment and overall
        % in the catch section.
        try
            plot(ax(ccc,cexp),ttx,...
                PSTHStruct.PSTHs(ccc).PSTHstack(1,:,cexp)./tpe,... ofst(cexp)+PSTHStruct.PSTHs(ccc).PSTHstack(1,:,cexp)./tpe,...
                'Color',cmap(1,:),'DisplayName',PSTHStruct.Trigger,...
                'LineStyle','-');
        catch
            plot(ax(ccc,cexp),ttx, auxPSTH(1,:)./tpe, 'Color',cmap(1,:),...
                'DisplayName',PSTHStruct.Trigger,'LineStyle','-');
        end
        ax(ccc,cexp).XAxis.Visible = 'off';
    end
    ax(ccc,cexp).XAxis.Visible = 'on';
end
fig_out = fig;
end