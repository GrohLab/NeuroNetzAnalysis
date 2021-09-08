function evolutionFigures = plotResponseEvolution(popResponse, varargin)
%PLOTRESPONSEEVOLUTION creates a plot with the given population mean
%response and the given parameters.
%   figure = plotResponseEvolution(popResponse, varargin)

%% Parsing inputs

p = inputParser;
% Required arguments
checkPopResponse = @(x) all([isstruct(x),...
    isfield(x,{'Mean','Confidence','FanoFactor','SEM'})]);
% Parameters
% Response/considered window duration
defRespWinDelta = 1; % seconds
checkRespWinDelta = @(x) all([isnumeric(x), x>0, numel(x) == 1]);
% Sampling frequency
%defFs = 1/3.33e-5;
%checkSampFreq = @(x) all([isnumeric(x), x>0, numel(x) == 1]);
% Trigger times
defTrigs = [];
checkTriggers = @(x) isnumeric(x) | isduration(x);
% Trial smoothing
defTrialNum = 1;

% Optional
defInMins = false;
checkInMins = @(x) all([isnumeric(x) | islogical(x), numel(x) == 1]);

% Error plotting
defError = 'Confidence';
checkError = @(x) any(strcmpi(x,{'Confidence','SEM'}));

p.addRequired('popResponse', checkPopResponse);
p.addParameter('windowDuration', defRespWinDelta, checkRespWinDelta);
%p.addParameter('samplingFrequency', defFs, checkSampFreq);
p.addParameter('triggerTimes', defTrigs, checkTriggers);
p.addParameter('trialBin', defTrialNum, @isPositiveIntegerValuedNumeric)
p.addOptional('inMinutes', defInMins, checkInMins);
p.addOptional('error', defError, checkError);

p.parse(popResponse, varargin{:})

popResponse = p.Results.popResponse;
%fs = p.Results.samplingFrequency;
winDelta = p.Results.windowDuration;
trigTms = p.Results.triggerTimes;
trialSmooth = p.Results.trialBin;
inMinFlag = p.Results.inMinutes;
errPlot = p.Results.error;
%% Validation for argument interactions
% Number of times per trigger 
Ntrigs = numel(trigTms);
% Number of triggers in modulations
NcondTrials = arrayfun(@(x) numel(x.Mean), popResponse);
NpopTrials = sum(NcondTrials);

if any(diff(NpopTrials) ~= 0)
    fprintf(2,'Different trial number per modulation?!\n');
    fprintf(2,'No way!\n');
    return
end
% If the number of triggers in the time axis is different than in
% modulations, exit.
if Ntrigs ~= 0 && any((NpopTrials - Ntrigs)~=0)
    fprintf(2,'The cardinality for the given trigger times do not match');
    fprintf(2,' the triggers in the population!\n')
    return
elseif Ntrigs == 0
    % No trial time axis given? Place the evolution as trial number.
    trigTms = 1:NpopTrials(1);
end

%%  Auxiliary variables
fnOpts = {'UniformOutput', false};
axOpts = {'NextPlot','add','Color','none','Box','off'};
expSubs = @(x) x(1):x(2);
[Ncond, Nmod] = size(popResponse); evolutionFigures = gobjects(Nmod,1);
muTrace = arrayfun(@(x) movmean(x.Mean, trialSmooth)./winDelta,...
    popResponse, fnOpts{:});
switch errPlot
    case {'Confidence','confidence'}
        confTrace =...
            arrayfun(@(x) movmean(x.Confidence, trialSmooth, 2)./winDelta,...
            popResponse, fnOpts{:});
    case {'SEM','sem'}
        confTrace =...
            arrayfun(@(x) movmean(x.Mean+[x.SEM;-x.SEM], trialSmooth, 2)./winDelta,...
            popResponse, fnOpts{:});
end
mxAx = cellfun(@(x) max(x(:)), confTrace);
mxAx = max(mxAx(:)); ax = gobjects(2,Nmod); mxAx = ceil(mxAx/10)*10;
histOpts = {'BinLimits', [0,mxAx],'BinWidth',2.5,...
    'Normalization','probability','EdgeColor','none',...
    'Orientation','horizontal'}; 
mxH = cellfun(@(x) max(histcounts(x, histOpts{1:6})), muTrace);
mxH = ceil(max(mxH(:))*1e2)/1e2;
trigTmsSubs = [[0;cumsum(NcondTrials(1:Ncond-1,1))]+1,cumsum(NcondTrials(:,1))];
xlabString = "Experimental time ";
if inMinFlag
    xlabString = xlabString + "[min]";
else
    xlabString = xlabString + "[s]";
end
strFrmt = 'Evolution for %s clusters'; 
clusterModName = ["potentiated";"depressed"];
permSubs = nchoosek(1:Ncond,2); Nperm = size(permSubs,1);
pM = arrayfun(@(x)...
    arrayfun(@(y) ranksum(muTrace{permSubs(y,:),x}), (1:Nperm)'),...
    1:Nmod, fnOpts{:});
while iscell(pM)
    pM = cat(2, pM{:});
end
mednMds = cellfun(@median, muTrace); signClrMap = [0.2*ones(1,3);0.8,0.17,0.05];
mednLnsX = mxH:mxH*0.08:mxH*(0.1*(Nperm -1)+1);
mxH = mednLnsX(end);
%% Plotting loop
for cmod = 1:Nmod
    evolutionFigures(cmod) = figure("Color",'w');
    % Evolution axis
    ax(1,cmod) = subplot(1,6,1:5,axOpts{:});
    % Mean trace for different conditions
    arrayfun(@(c)...
        plot(ax(1,cmod), trigTms(expSubs(trigTmsSubs(c,:))), muTrace{c,cmod},...
        'LineWidth',2), 1:Ncond);
    % Confidence range for the mean or SEM
    arrayfun(@(c)...
        plot(ax(1,cmod), trigTms(expSubs(trigTmsSubs(c,:))), confTrace{c,cmod}',...
        'LineStyle',':','Color',0.5*ones(1,3)), 1:Ncond);
    ylim(ax(1,cmod),[0, mxAx]); ylabel(ax(1,cmod), 'Response rate [Hz]')
    % Proportion axis
    yyaxis(ax(1,cmod), 'right'); arrayfun(@(c)...
        plot(ax(1,cmod), trigTms(expSubs(trigTmsSubs(c,:))), muTrace{c,cmod},...
        'LineStyle','none'), 1:Ncond); ylim(ax(1,cmod),[0, mxAx]);
    muAux = mean(muTrace{1, cmod});
    yticks(ax(1,cmod), muAux:muAux:mxAx); yticklabels(ax(1,cmod), 1:numel(yticks))
    set(get(ax(1,cmod),'YAxis'), 'Color',0.3*ones(1,3));
    if inMinFlag
        xticks(ax(1,cmod), 0:15*60:max(xlim(ax(1,cmod)))); 
        xticklabels(ax(1,cmod), xticks(ax(1,cmod))./60);
    end
    xlabel(ax(1,cmod), xlabString); 
    title(ax(1,cmod), sprintf(strFrmt, clusterModName(cmod)))
    % Histogram axis
    ax(2,cmod) = subplot(1,6,6,axOpts{:}); arrayfun(@(c)...
        histogram(ax(2,cmod), muTrace{c,cmod}, histOpts{:}), 1:Ncond);
    axis(ax(2,cmod), [0,mxH,0,mxAx]); set(get(ax(2,cmod),'YAxis'),...
        'Visible','off')
    % Plot significance on the y-distributions
    plot(ax(2,cmod),[zeros(1,Ncond);mxH*ones(1,Ncond)],...
        mednMds(:,cmod)'.*ones(2,Ncond),'LineStyle', '--',...
        'Color',0.55*ones(1,3))
    arrayfun(@(r) plot(ax(2,cmod), repmat(mednLnsX(r),2,1),...
        mednMds(permSubs(r,:),cmod),...
        'Color', signClrMap([1,2]*[pM(r,cmod)>0.05;pM(r,cmod)<0.05],:),...
        'Marker','+','UserData', pM(r,cmod),'LineWidth',1), 1:Nperm)
    xlabel(ax(2,cmod), 'Response proportion')
end

end