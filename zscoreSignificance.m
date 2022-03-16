function [signMat, evZmu, signMatpt, evZpt, alph] = ...
    zscoreSignificance(Counts, varargin)
%ZSCORESIGNIFICANCE computes the mean and standard deviation of the
%spontaneous window spiking activity to get the zscore of the spiking
%during the evoked window.
%   INPUTS
%       Counts - Cx2 cell array with C conditions and 2 time periods:
%                spontaneous (1) and evoked (2)
%   OPTINAL-INPUT name, value pair
%       'Alpha', alphaValue - a Ax1 vector of significance levels with 
%                             A >= 1 to compare the z-score values.
%   OUTPUTS 
%       signMat - Cx1 cell array containing a flag for z-score crossing
%                 different fixed significant levels.
%       evZmu - Cx1 cell array containing the mean z-score for all trials
%               per cluster
%       evZpt - Cx1 cell array containing all z-scores for each trial per
%               cluster
%       alph - Ax1 array contining the alpha values used
%% Parse input
% Must be cell, have two columns, and the number of clusters in all cells
% must be the same.
p = inputParser;
checkCounts = @(x) iscell(x) & size(x,2)==2 & ...
    ~std(cellfun(@(y) size(y,1), x), 0, "all");
% Alpha
defAlpha = 0.05; % 2.5% is 5% in zscore
checkAlpha = @(x) isnumeric(x) & isvector(x) & all(x>0 & x<0.2);
% Construction of input parser
p.addRequired('Counts', checkCounts);
p.addParameter('Alpha', defAlpha, checkAlpha);
% Parsing
p.parse(Counts, varargin{:});
% Saving results
Counts = p.Results.Counts;
alph = p.Results.Alpha;
%% Algorithm
% Auxiliary variables
fnOpts = {'UniformOutput', false}; mymean = @(x) mean(x,2,"omitnan");
thcmp = @(x,y) abs(x) > reshape(y,1,1,[]);
nDist = makedist('Normal', "mu", 0, "sigma", 1);
% Computing thresholds for user defined significance level
alph = 1 - alph;
signTh = arrayfun(@(z) fminbnd(@(y) ...
    norm(integral(@(x) nDist.pdf(x), -y, y) - z, 2), -3, 3), alph);
alph = 1 - alph;
% Mean and standard deviation of spontaneous counts
[~, muZ, sigZ] = cellfun(@(x) zscore(x, 1, 2), Counts(:,1), fnOpts{:});
% Z-score of evoked counts per trial (EVoked Zscore Per Trial)
evZpt = cellfun(@(x,y,z) (x - y)./z, Counts(:,2), muZ, sigZ, fnOpts{:});
signMatpt = cellfun(@(x) thcmp(x, signTh), evZpt, fnOpts{:});
% Trial proportion per significance level crossing per cluster per
% condition
signMat = cellfun(@(x) mean(x, 2, "omitnan"), signMatpt, fnOpts{:});
% Mean of the zscores for all trials
evZmu = cellfun(mymean, evZpt, fnOpts{:});
end