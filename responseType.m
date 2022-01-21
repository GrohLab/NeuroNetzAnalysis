function [mSubs, coeff] = responseType(PSTH, timeLapse, varargin)
%RESPONSETYPE computes the PCA on the given PSTHs and returns the groups
%which 
%   Detailed explanation goes here
%% Input parsing
% Required arguments: PSTH and timeLapse (timeLapse has no function as for
% now)
checkPSTH = @(x) all([isnumeric(x), ismatrix(x) || ndims(x) == 3]);
checkTimeLapse = @(x) all([isvector(x), numel(x) == 2, x(1) < x(2)]);
% Optional arguments: explained variance (explTh) default value 0.8
defExplTh = 4/5;
checkExplTh = @(x) all([isnumeric(x), x <= 1 & x > 0, isreal(x)]);

p = inputParser;
addRequired(p, 'PSTH', checkPSTH)
addRequired(p, 'timeLapse', checkTimeLapse)
addParameter(p, 'ExplainVarThreshold', defExplTh, checkExplTh)
parse(p, PSTH, timeLapse, varargin{:});
% Variable assignment
PSTH = p.Results.PSTH;
timeLapse = p.Results.timeLapse;
explTh = p.Results.ExplainVarThreshold;

%% Principal Component Factorization
[coeff, score, latent] = pca(PSTH'); 
explVar = cumsum(latent/sum(latent));
explVar2 = latent(2:end).\latent(1:end-1) - 1;
pcFlag = explVar >= explTh; pcFlag2 = explVar2 > 1;
Npc = find(~pcFlag(1:end-1) & pcFlag2, 1, "last");
Nv = find(pcFlag, 1, "first");

switch clustAlg
    case 'spect'
        mSubs = spectral_clustering(coeff, Nv, Npc, cpFlag);
    otherwise
        fprintf(1, "Unknown clustering algorithm!\n")
        return
end

end

function [mSubs, eVec] = spectral_clustering(X, k, Npc, plotFlag)
[mSubs, eVec] = spectralcluster(X(:,1:Npc), k,...
    'Distance', 'seuclidean');
if plotFlag
    figure;
    for ccl = unique(mSubs)'
        scatter3(eVec(mSubs == ccl,Npc-2), eVec(mSubs == ccl,Npc-1),...
            eVec(mSubs == ccl,Npc), "filled", "MarkerFaceAlpha", 0.4,...
            "DisplayName", num2str(ccl)); hold on;
    end
    text(eVec(:,Npc-2), eVec(:,Npc-1), eVec(:,Npc), pclID)
    lgnd = legend("show"); set(lgnd, "Box", "off", "Location", "best")
end
end