function [CV2, CVsqr, brstIdx, threshPrCl] =...
    getBurstInfo(spkSubs, fs, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Parsing input
p = inputParser;

checkFs = @(x) all([isnumeric(x), numel(x) == 1, x > 0]);

defNbin = 128;
checkNbin = @(x) all([isnumeric(x), numel(x) == 1, (round(x) - x) == 0,...
    x > 0]);

defMxDom = 1e3;
checkMxDom = @(x) all([isnumeric(x), 1/fs < x, x > 0, numel(x) == 1]);

addRequired(p, 'spkSubs', @iscell);
addRequired(p, 'fs', checkFs);
addOptional(p, 'Nbin', defNbin, checkNbin);
addOptional(p, 'MxDom', defMxDom, checkMxDom);


p.parse(spkSubs, fs, varargin{:});

spkSubs = p.Results.spkSubs;
fs = p.Results.fs;
Nbin = p.Results.Nbin;
mxDom = p.Results.MxDom;

%% Calculating ISIs and their log-distributions
fnOpts = {'UniformOutput', false};
%hOpts = {'Normalization', 'probability'};
Ncl = size(spkSubs,1);
% Nspks = cellfun(@numel, spkSubs);
% silentClusters = Nspks < 20;

% Preparing the logarithmic domain
% [isiCent, isiEdg, ~, logTau] = prepareLogBinEdges([1/fs, mxDom], Nbin);
[isiCent, ~, ~, logTau] = prepareLogBinEdges([1/fs, mxDom], Nbin);
% Getting the spike timing
spkTms = cellfun(@(x) x./fs, spkSubs, fnOpts{:});
% Nspks = cellfun(@numel, spkTms);
% Getting the inter-spike interval per cluster
isiTms = cellfun(@(x) diff(x), spkTms, fnOpts{:});
% Fitting a Gaussian Mixture Model to find clean peaks in the distribution
gmFt = cellfun(@(x) emforgmm(log10(x), 8), isiTms, fnOpts{:});
% Getting rid of overfitting artifacts
gmDelFlag = cellfun(@(x) (x(:,3)/logTau) < 1e-3, gmFt, fnOpts{:});
gmFt = cellfun(@(x,y) x(~y,:), gmFt, gmDelFlag, fnOpts{:});
% Computing the smooth distribution
px = cellfun(@(x) genP_x(x, isiCent), gmFt, fnOpts{:}); pxMat = cat(1, px{:});
% Extracting peaks and valleys (critical points)
[isiCp, isiSl] = getWaveformCriticalPoints(pxMat', 1/logTau);
% Adding offset to the found critical points and getting only the peaks
isiPks = cellfun(@(x,y) x(y<0) - log10(fs), isiCp(:,1), isiSl(:,1), fnOpts{:});
isiVls = cellfun(@(x,y) x(y>0) - log10(fs), isiCp(:,1), isiSl(:,1), fnOpts{:});
% Removing unrealistic peaks smaller than 1 ms
isiPks = cellfun(@(x) x(x>-3), isiPks, fnOpts{:});
isiVls = cellfun(@(x) x(x>(log10(1.5)-3)), isiVls, fnOpts{:});
isiPVl = arrayfun(@(x) interp1(isiCent, pxMat(x,:), isiPks{x}), (1:Ncl)', fnOpts{:});
isiPVz = cellfun(@(x,y) (x - mean(y))./std(y), isiPVl, px, fnOpts{:});
isiPks = cellfun(@(x,y) x(y>1), isiPks, isiPVz, fnOpts{:});
% isiPVl = arrayfun(@(x) interp1(isiCent, pxMat(x,:), isiPks{x}), (1:Ncl)', fnOpts{:});
% Interpolation of the peak value.

isiVVl = arrayfun(@(x) interp1(isiCent, pxMat(x,:), isiVls{x}), (1:Ncl)', fnOpts{:});
% isiVls = cellfun(@(x,y) x(y>0) - log10(fs), isiCp(:,1), isiSl(:,1), fnOpts{:});
% pts = cellfun(@(x,y) cat(2, x, y), isiPks, isiPVl, fnOpts{:});

retMps = cellfun(@(x) log10(cat(2, x(1:end-1), x(2:end))), isiTms, fnOpts{:});
Npks = cellfun(@numel, isiPks);

[CV2, CVsqr] = cellfun(@getCVsfromISIs, isiTms);

%%
% No easy way out... looping per cluster
valleyFlag = cellfun(@isempty, isiVls); Nvls = cellfun(@numel, isiVls);
threshPrCl = zeros(Ncl,1); brstIdx = threshPrCl;

    function thresh = getLowestValley(loc, val)
        [~, lwstVall] = min(val);
        thresh = evaluateAndAssignThreshVal(loc(lwstVall));
    end

    function thresh = evaluateAndAssignThreshVal(val)
        thresh = 0;
        if val <= -1
            thresh = val;
        end
    end

for ccl = 1:Ncl
    % If no valley was found, I'm assuming no bursts
    if ~valleyFlag(ccl)
        % If only one valley was found, select it as the threshold
        if Nvls(ccl) == 1 
            threshPrCl(ccl) = evaluateAndAssignThreshVal(isiVls{ccl});
        else
            % If more valleys were found, take a look at the peaks
            % If there was only one peak, take the lowest valley
            if Npks(ccl) == 1
                %[~, lowestValley] = min(isiVVl{ccl});
                %threshPrCl(ccl) = isiVls{ccl}(lowestValley);
                threshPrCl(ccl) = getLowestValley(isiVls{ccl}, isiVVl{ccl});
            else
                % Check the location of the valleys w.r.t. the peaks
                vallBetwnPksFlag = isiVls{ccl} < isiPks{ccl}';
                vall2right = sum(vallBetwnPksFlag);
                pks2right = sum(vallBetwnPksFlag,2);
                % If a valley to the left of the first peak is found, take
                % it as the threshold
                if vall2right(1) == 1
                    threshPrCl(ccl) = evaluateAndAssignThreshVal(isiVls{ccl}(1));
                elseif vall2right(1) > 1
                    bf1PkSub = pks2right == Npks(ccl);
                    threshPrCl(ccl) =...
                        getLowestValley(isiVls{ccl}(bf1PkSub),...
                        isiVVl{ccl}(bf1PkSub));
                else
                    if vall2right(2) == 1
                        threshPrCl(ccl) = evaluateAndAssignThreshVal(isiVls{ccl}(1));
                    else
                        bf1PkSub = pks2right == Npks(ccl) - 1;
                        threshPrCl(ccl) =...
                            getLowestValley(isiVls{ccl}(bf1PkSub),...
                            isiVVl{ccl}(bf1PkSub));
                    end
                end
            end
        end
    end
    if threshPrCl(ccl) < 0
        frstBrstSpk = retMps{ccl}(:,1) > threshPrCl(ccl) &...
            retMps{ccl}(:,2) < threshPrCl(ccl);
        tonicSpks = retMps{ccl}(:,1) > threshPrCl(ccl) &...
            retMps{ccl}(:,2) > threshPrCl(ccl);
        Nbrst = sum(frstBrstSpk); Ntonic = sum(tonicSpks);
        brstIdx(ccl) = Nbrst./(Nbrst + Ntonic);
    end
end


%%

% [matSubX, matSubY] = arrayfun(@(x) ind2sub([x,x], 1:x^2), Npks, fnOpts{:});
% matSubs = cellfun(@(x,y) cat(2, x(:), y(:)), matSubX, matSubY, fnOpts{:});
% retPts = cellfun(@(x,y) x(y), isiPks, matSubs, fnOpts{:});


% Difference of points
%dmCl = cellfun(@(x,y) distmatrix(x, y), retPts, retMps, fnOpts{:});


% Getting the distribution per cluster in log domain
% hisi = cellfun(@(x) histcounts(log10(x), isiEdg, hOpts{:}), isiTms, ...
%     fnOpts{:});
% hisi = cat(1,hisi{:});
% % Getting the peaks and valleys of the distribution bigger than 1 ms
% [hcp, xsl] = getWaveformCriticalPoints(hisi', 1/logTau);
% hcp = cellfun(@(x) x - log10(fs), hcp(:,1), fnOpts{:});
% msThFlags = cellfun(@(x) x > -3, hcp, fnOpts{:});
% hcp = cellfun(@(x,y) x(y), hcp, msThFlags, fnOpts{:});
% xsl = cellfun(@(x,y) x(y), xsl(:,1), msThFlags, fnOpts{:});
% hcp = cellfun(@(x,y) x(y<0), hcp, xsl, fnOpts{:});
% 
% ycp = arrayfun(@(x) interp1(isiCent, hisi(x,:), hcp{x}), (1:Ncl)',...
%     fnOpts{:});
% % Getting initial points for kmeans per cluster
% [~, ordSub] = cellfun(@(x) sort(x, 'descend'), ycp, fnOpts{:});
% 


end

