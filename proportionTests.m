function [p, h, fnc] = proportionTests(behRes, Na, varargin)

p = inputParser;
checkAlpha = @(x) isnumeric(x) & x > 0 & x < 1 & numel(x) == 1;

addRequired(p, 'behRes', @isstruct)
addRequired(p, 'Na', @isnumeric)
addParameter(p, 'alpha', 0.05, checkAlpha)

parse(p, behRes, Na, varargin{:})

behRes = p.Results.behRes;
Na = p.Results.Na;
alph = p.Results.alpha;

fnOpts = {'UniformOutput', false};
Nccond = numel(Na);
mt = arrayfun(@(c) [behRes(c).Results.MovProbability].*Na(c), 1:Nccond, ...
    fnOpts{:}); mt = cat(1, mt{:})';

    function ztest_prop()
        zDist = makedist("Normal", "mu", 0, "sigma", 1);
        p_hat = sum(mt,2)/sum(Na);
        Z = sum(mt./Na, 2) ./ sqrt(p_hat .* (1 - p_hat) * sum(1./Na));

        p = (1 - cdf(zDist, Z)) + cdf(zDist, -Z);
        h = p < alph;
    end

    function chiSqTest()
        mt_aux = mt';
        contingencyTbl = cat(3, mt_aux, Na' - mt_aux);
        contingencyTbl = permute(contingencyTbl, [1,3,2]);
        if any(contingencyTbl<=5, "all") && ...
                all(size(contingencyTbl,[1,2])==[2,2])
            [h, p] = arrayfun(@(ct) fishertest(round(contingencyTbl(:,:,ct))), ...
                1:size(contingencyTbl,3));
        else
            sumCol = sum(contingencyTbl, 1); sumCol = squeeze(sumCol);
            expectedTrials = arrayfun(@(r) (sumCol(:,r) * Na)./ ...
                sum(contingencyTbl(:,:,r), "all"), 1:size(sumCol,2), fnOpts{:});
            expectedTrials = cat(3, expectedTrials{:});
            Chi2 = ((contingencyTbl - expectedTrials).^2)./expectedTrials;
            Chi2 = sum(Chi2, [1,2]);
            p = chi2pdf(squeeze(Chi2), prod(size(contingencyTbl, [1,2])-1));
            h = p < alph;
        end
    end
        
end