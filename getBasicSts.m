function [stStruct, spkDom] = getBasicSts(spkCell, cStruct, Na, res)
%GETBASICSTS computes basic statistics on the spike trains for all trials
%in the relative spike times structure. 
%   Detailed explanation goes here
%%
if ~exist("res", "var") || ~(isnumeric(res) && numel(res) == 1 && res > 0)
    res = 5e-4;
end
foSts = @(x) [median(x), mean(x), std(x)];
expnd = @(x, r) linspace(x(1), x(end), 1+ceil(diff(x([1,end]))/res));
fnOpts = {'UniformOutput', false};
erOpts = {'ErrorHandler', @emptySpkTrain};
ksOpts = {'Function','pdf'};
%% 

spkDom = expnd(cStruct.Response_window_s, res); N = length(spkDom);
% Median, Mean, and STD
sFOS = cellfun(@(x) cellfun(foSts, x, fnOpts{:}), spkCell, fnOpts{:});
sFOS = cellfun(@(x) cat(1, x{:}), sFOS, fnOpts{:});
% Empty spike trains
eptyFlag = cellfun(@(x) cellfun(@(y) isempty(y), x), spkCell, fnOpts{:});
% Inverse cumulative density function of the spikes
sH = cellfun(@(x) cellfun(@(y) ksdensity(y, spkDom, ksOpts{:}), x, ...
    fnOpts{:}, erOpts{:}), spkCell, fnOpts{:});
for ca = 1:length(sH)
    for cb = reshape(find(eptyFlag{ca}),1,[])
        sH{ca}{cb} = repmat(sH{ca}{cb}, 1, N);
    end
end
% Quartile positions (0.25, 0.5, 0.75)
Qv = cellfun(@(a) arrayfun(@(z) cellfun(@(y) fzero(@(x) ...
    interp1(spkDom, cumsum(y./sum(y))-z, x), spkDom([1,end])), a, ...
    erOpts{:}), 0.25:0.25:0.75, fnOpts{:}), sH, fnOpts{:});
Qv = cellfun(@(x) cat(2, x{:}), Qv, fnOpts{:});
% Mahalanobis distance in matrix form
sMh = cellfun(@(x) cellfun(@(y) squareform(pdist(y(:), "mahalanobis")), ...
    x, fnOpts{:}, erOpts{1},@covError), spkCell, fnOpts{:});
% Saving results
stStruct = struct('FOstats', sFOS, 'PDF', sH,'Quartiles', Qv, ...
    'MahalDist', sMh);
end

function fsH = emptySpkTrain(S, varargin)
fsH = 0;
end

function sMh = covError(S, varargin)
sMh = [];
end