function ax = plotPopPSTH(PSTHStruct, varargin)
%PLOTPOPPSTH This function accepts the PSTH structure created by getPopPSTH
%and plots all PSTHs of all experiments together with a population PSTH at
%the bottom of the figure. It can plot the PSTH structure in a given axis
%or create a new figure to plot this results.
%   Emilio Isaías-Camacho @ GrohLab 2019

%% Input validation
if isstruct(PSTHStruct) && isfield(PSTHStruct,'TimeAxis')
    p = inputParser;
    % Bin size validation
    defaultBin = PSTHStruct.BinSize;
    small_limit = mean(diff(PSTHStruct.TimeAxis));
    attributes = {'scalar','>',small_limit,...
        '<=', sum(abs([PSTHStruct.TimeAxis(1),PSTHStruct.TimeAxis(end)]))*0.1};
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
    % Actual validation
    p.parse(PSTHStruct, varargin{:})
else
    fprintf('The PSTH structure is either corrupt or hasn''t been created\n')
    return
end
%% Initialization
% Bin size
binsize = p.Results.BinSize;
ttx = PSTHStruct.TimeAxis;
btx = ttx(1):binsize:ttx(end);
% Figure
fig = p.Results.Figure;
fig.Children = axes();
% Number of conditions combinations
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

%% Plotting procedure
[Nv, Ns, Nx] = size(PSTHStruct.PSTHs(1).PSTHstack);
ax = gobjects(Ncc,Nx);
fs = 1/mean(diff(ttx));
cmap = winter(Ncc-1);
cmap = [0.1,0.1,0.1;cmap];
fig.Visible = 'on';
% ofst = (0:Nx-1) * 1.2;
Nb = numel(btx);
for ccc = 1:Ncc
    % Computing the necessary counters and subscripts to manage the
    % multiple subplots in the graphs.
    ccol = mod(ccc-1,Ncol)+1;
    for cexp = 1:Nx
        crow = ((ccc - ccol)/Ncol);
        crowFig = (crow*(Nx)) + cexp;
        sIdx = sub2ind([Nrow*(Nx),Ncol],crowFig,ccol);
        ax(ccc,cexp) = subplot(Nrow*(Nx),Ncol,sIdx);
        ax(ccc,cexp).Box = 'off';
        ax(ccc,cexp).XAxis.Visible = 'off';
        ax(ccc,cexp).Parent = fig;
        hold(ax(ccc,cexp),'on')
        % Binning the signals using the _binsize_ variable
        [bPSTH, binWidth] =...
            binTime(PSTHStruct.PSTHs(ccc).PSTHstack(2:Nv,:,cexp),binsize,fs);
        
        % Removing those variables which didn't occure at all
        emptyFlag = sum(bPSTH,1) == 0;
        bPSTH = bPSTH(:,~emptyFlag);
        nNv = size(bPSTH,2);
        lbs = PSTHStruct.SignalIDs([false,~emptyFlag]);
        cSubs = find([false,~emptyFlag]);
        % Normalization values
        tpe = PSTHStruct.PSTHs(ccc).TrialsPerExperiment(cexp);
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
                    'DisplayName',lbs{cvar})
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
        plot(ax(ccc,cexp),ttx,...
            PSTHStruct.PSTHs(ccc).PSTHstack(1,:,cexp)./tpe,... ofst(cexp)+PSTHStruct.PSTHs(ccc).PSTHstack(1,:,cexp)./tpe,...
            'Color',cmap(1,:),'DisplayName',PSTHStruct.Trigger,...
            'LineStyle','-')
    end
    ax(ccc,cexp).XAxis.Visible = 'on';
end
end

% p = inputParser;
% charchk = {'char'};
% numchk = {'numeric'};
% nempty = {'nonempty'};
%
% addRequired(p,'shape',@(x)validateattributes(x,charchk,nempty))
% addRequired(p,'dim1',@(x)validateattributes(x,numchk,nempty))
% addOptional(p,'dim2',1,@(x)validateattributes(x,numchk,nempty))
% parse(p,shape,dim1,varargin{:})
