function [p, h] = proportionTests(behRes, Na, varargin)

ip = inputParser;
checkAlpha = @(x) isnumeric(x) & x > 0 & x < 1 & numel(x) == 1;

defTest = "ztest";
testPoss = [defTest, "chi"];
istxt = @(x) ischar(x) | isstring(x);
checkTest = @(x) istxt(x) & any(strcmpi(x, testPoss));

addRequired(ip, 'behRes', @isstruct)
addRequired(ip, 'Na', @isnumeric)
addParameter(ip, 'alpha', 0.05, checkAlpha)
addParameter(ip, 'test', defTest, checkTest)

parse(ip, behRes, Na, varargin{:})

behRes = ip.Results.behRes;
Na = ip.Results.Na;
alph = ip.Results.alpha;

fnOpts = {'UniformOutput', false};
Nccond = numel(Na);
mt = arrayfun(@(c) [behRes(c).Results.MovProbability].*Na(c), 1:Nccond, ...
    fnOpts{:}); mt = cat(1, mt{:})';



    function ztest_prop()
        zDist = makedist("Normal", "mu", 0, "sigma", 1);
        p_hat = sum(mt,2)/sum(Na);
        Z = sum(mt./Na, 2) ./ sqrt(p_hat .* (1 - p_hat) * sum(1./Na));
        % Tow-tail test
        p = (1 - cdf(zDist, Z)) + cdf(zDist, -Z);
        h = p < alph;
    end

    function chiSqTest()
        mt_aux = mt';
        contingencyTbl = cat(3, mt_aux, Na' - mt_aux);
        contingencyTbl = permute(contingencyTbl, [1,3,2]);
        [r,c] = size(contingencyTbl, [1,2]);
        if any(contingencyTbl<=5, "all") && ...
                all([r,c]==2, "all")
            [h, p] = arrayfun(@(ct) fishertest(round(contingencyTbl(:,:,ct))), ...
                1:size(contingencyTbl,3));
        else
            sumCol = sum(contingencyTbl, 1); %sumCol = squeeze(sumCol);
            expectedTrials = arrayfun(@(r) (Na' * squeeze(sumCol(:,:,r)))./ ...
                sum(contingencyTbl(:,:,r), "all"), 1:size(sumCol,3), fnOpts{:});
            expectedTrials = cat(3, expectedTrials{:});
            Chi2 = ((contingencyTbl - expectedTrials).^2)./expectedTrials;
            Chi2 = sum(Chi2, [1,2]);
            p = chi2pdf(squeeze(Chi2), prod(size(contingencyTbl, [1,2])-1));
            h = p < alph;
        end
    end
        
end