function [mSubs, coeff, Npc, Nv] = responseType(PSTH, timeLapse, varargin)
%RESPONSETYPE computes the PCA on the given PSTHs and returns the groups
%which 
%   Detailed explanation goes here

%% Input parsing
% Required arguments: PSTH and timeLapse (timeLapse has no function as for
% now)
checkPSTH = @(x) all([isnumeric(x), ismatrix(x)]);
checkTimeLapse = @(x) all([isvector(x), numel(x) == 2, x(1) < x(2)]);

% Optional arguments: explained variance (explTh) default value 0.8
defExplTh = 4/5;
checkExplTh = @(x) all([isnumeric(x) & isscalar(x), x <= 1 & x > 0]);

% Algorithms: spect, heirarchical, kmeans
defAlg = 'spect';
validAlg = {'spect', 'hierarchical', 'kmeans'};
checkAlg = @(x) any(validatestring(x, validAlg));

p = inputParser;
addRequired(p, 'PSTH', checkPSTH)
addRequired(p, 'timeLapse', checkTimeLapse)
addParameter(p, 'ExplainVarThreshold', defExplTh, checkExplTh)
addParameter(p, 'Algorithm', defAlg, checkAlg)
addParameter(p, 'Plot', false, @(x) islogical(x) && numel(x) == 1)

parse(p, PSTH, timeLapse, varargin{:});
% Variable assignment
PSTH = p.Results.PSTH;
timeLapse = p.Results.timeLapse;
explTh = p.Results.ExplainVarThreshold;
clustAlg = p.Results.Algorithm;
cpFlag = p.Results.Plot;

%% Principal Component Factorization
[Ncl, Nts] = size(PSTH);
[coeff, score, latent] = pca(PSTH'); 
explVar = cumsum(latent/sum(latent));
explVar2 = latent(2:end).\latent(1:end-1) - 1;
pcFlag = explVar >= explTh; pcFlag2 = explVar2 > 1;
Npc = find(~pcFlag(1:end-1) & pcFlag2, 1, "last");
Nv = find(pcFlag, 1, "first");

mSubs = zeros(Ncl,1);
%% Clustering algorithm
switch clustAlg
    case 'spect'
        mSubs = spectral_clustering(coeff, Nv, Npc, cpFlag);
    case 'heirarchical'
        mSubs = heirarchical_clustering(coeff, Nv, Npc, cpFlag);
    case 'kmeans'
        [mSubs, mC] = kmeans_clustering(coeff, Nv, Npc, cpFlag);
    case 'kmedoids'
        [mSubs, mC] = kmedoids_clustering(coeff, Nv, Npc, cpFlag);
    otherwise
        fprintf(1, "Unknown clustering algorithm!\n")
        return
end

end
%% Local functions
function [mSubs, eVec] = spectral_clustering(X, k, Npc, plotFlag)
[mSubs, eVec, eVal] = spectralcluster(X(:,1:Npc), k,...
    'Distance', 'seuclidean');
if plotFlag
    figure;
    for ccl = unique(mSubs)'
        scatter3(eVec(mSubs == ccl,Npc-2), eVec(mSubs == ccl,Npc-1),...
            eVec(mSubs == ccl,Npc), "filled", "MarkerFaceAlpha", 0.4,...
            "DisplayName", num2str(ccl)); hold on;
        mC = mean(eVec(mSubs == ccl,:), 'omitnan');
        text(mC(:,Npc-2), mC(:,Npc-1), mC(:,Npc), num2str(ccl))
    end
    
    lgnd = legend("show"); set(lgnd, "Box", "off", "Location", "best")
end
end

function mSubs = heirarchical_clustering(X, k, Npc, plotFlag)
dx = pdist(X(:,1:k), "mahalanobis");
lx = linkage(dx);
mSubs = cluster(lx, 'maxclust', Npc);
if plotFlag
    figure; dendrogram(lx, max(mSubs));
end
end

function [mSubs, mC] = kmeans_clustering(X, k, Npc, plotFlag)
[mSubs, mC] = kmeans(X(:,1:Npc), k, "Distance", "sqeuclidean");
if plotFlag
    figure;
clrMap = lines(max(mSubs));
scObj = gobjects(max(mSubs),1);
for ccl = unique(mSubs)'
    scatter3(coeff(mSubs == ccl, Npc-2), coeff(mSubs == ccl, Npc-1),...
        coeff(mSubs == ccl, Npc), "filled", "MarkerFaceAlpha", 0.4,...
        "MarkerFaceColor", clrMap(ccl,:), "MarkerEdgeColor", "none",...
        "DisplayName", num2str(ccl));
    hold on;
    scObj(ccl) = scatter3(mC(ccl,1), mC(ccl,2), mC(ccl,3), "filled", ...
        "DisplayName", ['Centre ' num2str(ccl)],...
        "MarkerFaceColor", clrMap(ccl,:));
end
lgnd = legend(scObj); set(lgnd, "Box", "off", "Location", "best")
end

end

function [mSubs, mC] = kmedoids_clustering(X, k, Npc, plotFlag)
[mSubs, mC] = kmedoids(X(:,1:Npc), k, "Distance", "sqeuclidean");
if plotFlag
    figure;
clrMap = lines(max(mSubs));
scObj = gobjects(max(mSubs),1);
for ccl = unique(mSubs)'
    scatter3(coeff(mSubs == ccl, Npc-2), coeff(mSubs == ccl, Npc-1),...
        coeff(mSubs == ccl, Npc), "filled", "MarkerFaceAlpha", 0.4,...
        "MarkerFaceColor", clrMap(ccl,:), "MarkerEdgeColor", "none",...
        "DisplayName", num2str(ccl));
    hold on;
    scObj(ccl) = scatter3(mC(ccl,1), mC(ccl,2), mC(ccl,3), "filled", ...
        "DisplayName", ['Centre ' num2str(ccl)],...
        "MarkerFaceColor", clrMap(ccl,:));
end
lgnd = legend(scObj); set(lgnd, "Box", "off", "Location", "best")
end

end