function [corrMat, MiCorr] = getCorrInConditions(spkSubs,condArray,fs,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Input parsing
p = inputParser;

checkCondArray = @(x) all([isstruct(x), all(isfield(x, {'name','Triggers'}))]);
checkFS = @(x) all([isnumeric(x), numel(x) == 1, x>0]);

defBinWidth = 1.5; % Seconds
checkBinWidth = @(x) all([isnumeric(x), x > 0, numel(x) == 1]);

p.addRequired('spkSubs',@iscell);
p.addRequired('condArray',checkCondArray);
p.addRequired('fs',checkFS);
p.addParameter('binWidth', defBinWidth, checkBinWidth);

p.parse(spkSubs, condArray, fs, varargin{:})

spkSubs = p.Results.spkSubs;
condArray = p.Results.condArray;
fs = p.Results.fs;
bnSz = p.Results.binWidth;

%% Corrlation matrices per condition duration in block (maybe as trial later)
fnOpts = {'UniformOutput', 0}; Ncond = numel(condArray);
histOpts = {'Normalization', 'countdensity', 'BinWidth', bnSz, 'BinLimits'};
getMI = @(x, d) diff(x,1,d)./sum(x,d,'omitnan');

% Getting the time limits for the conditions in block
condLims = arrayfun(@(x) x.Triggers([1,end])./fs, condArray', fnOpts{:});
condLims = cat(1,condLims{:}).*[0.95, 1.05];

bnFr = arrayfun(@(y)...
    cellfun(@(x) histcounts(x./fs, histOpts{:}, condLims(y,:)),...
    spkSubs, fnOpts{:}), (1:Ncond)', fnOpts{:});
bnFr = cellfun(@(x) cat(1, x{:}), bnFr, fnOpts{:});
corrMat = cellfun(@(x) corrcoef(x'), bnFr, fnOpts{:});
corrMat = cat(3, corrMat{:}); condPerm = nchoosek(1:Ncond,2); 
Nperm = size(condPerm,1);

MiCorr = arrayfun(@(x) getMI(corrMat(:,:,condPerm(x,:)),3), (1:Nperm)', fnOpts{:}); 

end

% Code to compare spontaneous firing in two different conditions

% sponFr = arrayfun(@(x) arrayfun(@(z)  histcounts(condSpks{x,z},...
% 'BinLimits', condLims(z,:), 'BinWidth', 1, 'Normalization',...
% 'countdensity'), (1:Ncond) ,'UniformOutput', 0),
% (1:size(condSpks,1))',... 
% 'UniformOutput', 0);

% sponFr = cat(1, sponFr{:})