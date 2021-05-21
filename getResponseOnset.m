function [respOnset, qVals, fspkPcond] = getResponseOnset(relSpkTms, respTmWin,...
    varargin)
%GETRESPONSEONSET 
%   Detailed explanation goes here
%Emilio Isaias-Camacho @GrohLab 2021
%% Parsing inputs
p = inputParser;
% Required arguments
checkSpkTms = @(x) all([isstruct(x), isfield(x, {'name','SpikeTimes'})]);
checkTmWin = @(x) all([numel(x) == 2, x(1) < x(2)]);
% Parameters
% Bin size
defBnSz = 1e-3;
checkBinSize = @(x) all([isnumeric(x), x > 0, numel(x) == 1]);
% Threshold for response distribution
defAlpha = 0.05;
checkAlpha = checkBinSize;

p.addRequired('relSpkTms',checkSpkTms);
p.addRequired('respTmWin',checkTmWin);
p.addParameter('binWidth',defBnSz, checkBinSize);
p.addParameter('alpha',defAlpha, checkAlpha);

p.parse(relSpkTms, respTmWin, varargin{:});

relSpkTms = p.Results.relSpkTms;
tmWin = p.Results.respTmWin;
bnSz = p.Results.binWidth;
alph = p.Results.alpha;

%% Auxiliary variables
%{
     function [dataTmWin] = getTimeLimitsFromData(relSpkTms)
        auxIdx = 1;
        dataTmWin = zeros(2,1);
        for func = {@min, @max}
            auxTm = arrayfun(@(x) cellfun(func{:},...
                relSpkTms(x).SpikeTimes, fnOpts{:}), (1:Ncond)', fnOpts{:});
            while iscell(auxTm)
                auxTm = cat(2, auxTm{:});
            end
            dataTmWin(auxIdx) = func{:}(auxTm); auxIdx = auxIdx + 1;
        end
    end
%}
% Cell- and arrayfun packing results into a cellarray
fnOpts = {'UniformOutput', false};
% PSTH per cluster normalized as proportion
psthOpts = {'BinLimits',tmWin,'BinWidth',bnSz};
% Bin width for calculating the weighted average first spike.
frstSpkBnSz = 7.5e-4;
% Histogam options for the first spike histogram
fskOpts = {'Normalization','probability','BinWidth', frstSpkBnSz}; % 0.75 ms
% Number of clusters
Ncl = size(relSpkTms(1).SpikeTimes,1);

%% Getting the response onset per cluster
% Concatenating all spikes from all trials per cluster per condition
spkClAllTrials = arrayfun(@(x)...
    arrayfun(@(y) cat(2, x.SpikeTimes{y,:}), (1:Ncl)', fnOpts{:}),...
    relSpkTms, fnOpts{:});
% Computing PSTH per cluster per condition
[clPsth, bnEdg] = cellfun(@(x)...
    cellfun(@(y) histcounts(y, psthOpts{:}) ,x ,fnOpts{:}),...
    spkClAllTrials, fnOpts{:}); bnEdg = bnEdg{1,1}{1};
% Rearranging the PSTHs into a matrix
clPsth = cellfun(@(x) cat(1, x{:}), clPsth, fnOpts{:});
clPsth = cellfun(@(x) x./sum(x,2), clPsth, fnOpts{:});
cumPsth = cellfun(@(x) cumsum(x,2,'omitnan'), clPsth, fnOpts{:});
% Computing the time axis
deltaBin = mean(diff(bnEdg)); centEdg = bnEdg([1,end]) + [1,-1]*(deltaBin/2);
psthCent = centEdg(1):deltaBin:centEdg(2);
% Crossing of the response with alpha value
alphaCross = cellfun(@(y)...
    arrayfun(@(x) find(y(x,:) < alph, 1, 'last'), (1:Ncl)', fnOpts{:}),...
    cumPsth, fnOpts{:});
% Adjusting unexistant crossings with the sum (empty crossings) and with
% the anonymous function (all signal is below the threshold).
alphaCross = cellfun(@(x)...
    cellfun(@(y) sum([isempty(y),y],'omitnan'),x , fnOpts{:}),...
    alphaCross, fnOpts{:});
adjustSub = @(x,y) repmat(([1, 1;-1, 0]*[x==numel(y);x<numel(y)])'*[x;1],2,1)+[0;1];
% Linear model for calculating the x value for y = alpha
intrpLine = cellfun(@(x,y)...
    arrayfun(@(u)...
    fit_poly(psthCent(adjustSub(x{u},psthCent)),...
    y(u,adjustSub(x{u},psthCent)), 1), (1:Ncl)', fnOpts{:}),...
    alphaCross, cumPsth, fnOpts{:});
% Readjusting the linear equations to remove one cellfun.
intrpLine = cellfun(@(x) cat(2, x{:}), intrpLine, fnOpts{:});
respOnset = cellfun(@(x) (alph - x(2,:))./x(1,:), intrpLine, fnOpts{:});
% Getting default alpha threshold crossing.
[~,~,qVals] = cellfun(@(x) exponentialSpread(x, psthCent, tmWin),...
    clPsth, fnOpts{:});

%% Getting the first spike 
% After the stimulus onset per cluster per trial per condition.
clFspk = arrayfun(@(y)...
    cellfun(@(x) x(find(x>0,1,'first')), y.SpikeTimes, fnOpts{:}),...
    relSpkTms, fnOpts{:});
% Rearranging all first spikes in a vector per cluster per condition.
clFspk = cellfun(@(x)...
    arrayfun(@(y) cat(2,x{y,:}), (1:Ncl)', fnOpts{:}),...
    clFspk, fnOpts{:});
% Creating the distribution per cluster per condition
[fspkHist, bnEdg] = cellfun(@(x)...
    cellfun(@(y) histcounts(y, fskOpts{:}), x, fnOpts{:}),...
    clFspk, fnOpts{:});
bnCent = cellfun(@(x)...
    cellfun(@(y) (y(1)+(frstSpkBnSz/2)):frstSpkBnSz:(y(end)-(frstSpkBnSz/2)),...
    x, fnOpts{:}),bnEdg, fnOpts{:});
% Computing the dot product for weighted average for the first spike
fspkPcond = cellfun(@(x,y)...
    cellfun(@(u,v) u*v', x, y, fnOpts{:}),...
    fspkHist, bnCent, fnOpts{:});
% Rearranging the results
fspkPcond = cellfun(@(x) cat(1,x{:}), fspkPcond, fnOpts{:});
end

