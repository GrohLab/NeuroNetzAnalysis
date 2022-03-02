function [signMat, evZ, evZpt, alph] = zscoreSignificance(Counts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fnOpts = {'UniformOutput', false}; mymean = @(x) mean(x,2,"omitnan");
thcmp = @(x,y) x > reshape(y,1,1,[]);
nDist = makedist('Normal', "mu", 0, "sigma", 1);
% Replace for user input thresholds
alph = [0.95, 0.98, 0.99, 0.995, 0.998, 0.999]; 
signTh = arrayfun(@(y) fminbnd(@(x) norm(nDist.cdf(x) - y, 2), 1, 5),...
    alph);
[~, muZ, sigZ] = cellfun(@(x) zscore(x, 1, 2), Counts(:,1), fnOpts{:});
evZpt = cellfun(@(x,y,z) (x - y)./z, Counts(:,2), muZ, sigZ, fnOpts{:});
evZ = cellfun(mymean, evZpt, fnOpts{:});
signMat = cellfun(@(x) thcmp(x, signTh), evZ, fnOpts{:});
end