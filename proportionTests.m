function [p, h, fnc] = proportionTests(behRes, Na, varargin)

p = inputParser;
checkAlpha = @(x) isnumeric(x) & x > 0 & x < 1 & numel(x) == 1;

addRequired(p, 'behRes', @isstruct)
addRequired(p, 'Na', @isnumeric)
addParameter(p, 'alpha', 0.05, checkAlpha)

parse(p, behRes, Na, varargin{:})

behRes = p.Results.behRes;
Na = p.Results.Na;
alpha = p.Results.alpha;

fnOpts = {'UniformOutput', false};
Nccond = numel(Na);
mt = arrayfun(@(c) [behRes(c).Results.MovProbability].*Na(c), 1:Nccond, ...
            fnOpts{:}); mt = cat(1, mt{:})';
    function ztest_prop()
        zDist = makedist("Normal", "mu", 0, "sigma", 1);
        p_hat = sum(mt,2)/sum(Na);
        Z = sum(mt./Na, 2) ./ sqrt(p_hat .* (1 - p_hat) * sum(1./Na));

        p = (1 - cdf(zDist, Z)) + cdf(zDist, -Z);
        h = p < alpha;
    end
    function chi2test()
        mt_aux = mt';
        contingencyTbl = cat(3, mt_aux, Na' - mt_aux);
        sumCol = sum(contingencyTbl, 1); sumCol = squeeze(sumCol);
        expectedTrials = arrayfun(@(r) (sumRow(:,r) * Na)./sum(Na), 1:size(sumRow,2), fnOpts{:});
        contingencyTbl = permute(contingencyTbl, [1,3,2]);
        
        mt = arrayfun(@(c) [behRes(c).Results.MovProbability].*Na(c), 1:Nccond, ...
fnOpts{:}); mt = cat(1, mt{:})';
mt
mt_aux = mt';
mt_aux
Na' - mt_aux
contingencyTbl = cat(3, mt_aux, Na' - mt_aux);
contingencyTbl
sum(contingencyTbl, 1)
sum(contingencyTbl, 2)
sum(contingencyTbl, 3)
sumRow = sum(contingencyTbl, 3)
sumCol = sum(contingencyTbl, 1);
sumCol
sum(contingencyTbl, [1,2])
sum(contingencyTbl, [1,3])
size(sumRow
size(sumRow)
size(sumCol)
size(contingencyTbl)
sumRow
sumCol * sumRow
pagetimes(sumCol, sumRow)
pagemtimes(sumCol, sumRow)
sumCol
sumCol = squeeze(sumCol);
sumCol
size(sumCol)
sumRow
mt_aux
contingencyTbl
mt_aux
sumRow(:,1) * Na
(sumRow(:,1) * Na)./sum(Na)
(sumRow(:,1) * reshape(Na, 1,1,size(Na,2))./sum(Na)
(sumRow(:,1) * reshape(Na, 1,1,size(Na,2))./sum(Na))
(sumRow * reshape(Na, 1,1,size(Na,2))./sum(Na))
tensorprod(sumRow, Na)
arrayfun(@(r) (sumRow(:,r) * Na)./sum(Na), 1:size(sumRow,2), fnOpts{:})
expectedTrials = arrayfun(@(r) (sumRow(:,r) * Na)./sum(Na), 1:size(sumRow,2), fnOpts{:});
contingencyTbl
transp(contingencyTbl)
pagetransp(contingencyTbl)
permute(contingencyTbl, [1,3,2])
contingencyTbl
contingencyTbl = permute(contingencyTbl, [1,3,2]);
contingencyTbl
expectedTrials
expectedTrials = cat(3, expectedTrials{:});
expectedTrials
resChi = ((contingencyTbl - expectedTrials).^2) ./ expectedTrials
sum(resChi, [1,2])
chi2pdf(sum(resChi, [1,2]), 1)
chi2pdf(sum(resChi, [1,2]), 1) < 0.05
    
    end
end