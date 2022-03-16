function [rclIdx, medH, zH] = getSignificantFlags(Results, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% input parsing
p = inputParser;
checkRes = @(x) isstruct(x) & all(contains(fieldnames(x), ...
    {'Combination', 'Activity'}));
% Alpha
defAlpha = 0.05; % 2.5% is 5% in zscore
checkAlpha = @(x) isnumeric(x) & numel(x)==1 & all(x>0 & x<0.2);
% Z-proportions
defZProp = 0.5;
checkZProp =  @(x) isnumeric(x) & numel(x)==1 & all(x>0.5 & x<=1);
% Construction of input parser
p.addRequired('Results', checkRes)
p.addParameter('Alpha', defAlpha, checkAlpha);
p.addParameter('ZProp', defZProp, checkZProp)
% Parsing
p.parse(Results, varargin{:});
Results = p.Results.Results;
alph = p.Results.Alpha;
zProp = p.Results.ZProp;
%% Estimation of number of conditions
fnOpts = {'UniformOutput', false};
Nres = size(Results, 2); auxM = triu(squareform(1:Nres, "tomatrix"));
auxM(:,1) = []; auxM(end,:) = []; Nccond = size(auxM,1); 
indCondSubs = auxM(:,Nccond);
% Median test: is spontaneous spike count median different from the median
% in evoked?
medH = cell2mat(cellfun(@(y) y(1).Pvalues,...
    arrayfun(@(x) x.Activity, Results(indCondSubs), fnOpts{:}),...
    fnOpts{:})) < alph;
% Z-score test: Proportion of trials passing significance level
zH = cell2mat(cellfun(@(y) y(1).ZSignificance,...
    arrayfun(@(x) x.Activity, Results(indCondSubs), fnOpts{:}),...
    fnOpts{:})) > zProp;
%% Passing both tests: 
% Different median between spontaneous and evoked, and z-score higher than
% spontaneous for, at least, more than half trials
rclIdx = medH & zH;
end
