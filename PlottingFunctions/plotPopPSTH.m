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
ax = gobjects(Ncc,1);
[Nv, Ns, Nx] = size(PSTHStruct.PSTHs(1).PSTHstack);
fs = 1/mean(diff(ttx));
cmap = jet(Ncc);
fig.Visible = 'on';
for csp = 1:Ncc
    ax(csp) = subplot(Nrow,Ncol,csp,fig.Children);
    hold(ax(csp),'on')
    for cexp = Nx:-1:1
        % Binning the signals using the _binsize_ variable
        [bPSTH, binWidth] =...
            binTime(PSTHStruct.PSTHs(csp).PSTHstack(2:Nv,:,cexp),binsize,fs);
        % Removing those variables which didn't occure at all
        emptyFlag = sum(bPSTH,1) == 0;
        bPSTH = bPSTH(:,~emptyFlag);
        nNv = size(bPSTH,2);
        lbs = PSTHStruct.SignalIDs([false,~emptyFlag]);
        cSubs = find([false,~emptyFlag]);
        for cv = 1:nNv
            bar(ax(csp),btx,(cexp+1)*bPSTH(:,cv),'EdgeColor','none',...
                'FaceAlpha',0.3,'FaceColor',cmap(cSubs(cv),:),...
                'DisplayName',lbs{cv},'BaseValue',cexp)
        end
        plot(ax(csp),ttx,PSTHStruct.PSTHs(csp).PSTHstack(1,:,cexp),...
            'Color',cmap(1,:),'DisplayName',PSTHStruct.Trigger)
        legend show
    end
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
