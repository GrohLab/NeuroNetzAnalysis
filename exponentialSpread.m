function [mdls, r2, qVals, qDiff] =...
    exponentialSpread(signMat, tx, tmWin)
%EXPONENTIALSPREAD measures where 50% of the data is sitting in a
%semi-poisson distribution and the exponential decay of the 1 - cumulative
%sum of the given signals.
%   Detailed explanation goes here
Ns = size(signMat,1);
wIdx = tx >= tmWin(1) & tx <= tmWin(2);
shTx = reshape(tx(wIdx), sum(wIdx), 1);
wSum = sum(signMat(:,wIdx),2); zSum = wSum ~= 0;
if ~all(zSum)
    warning(['Some clusters have no spikes between ',num2str(tmWin(1)*1e3),...
        ' and ', num2str(tmWin(2)*1e3),' ms relative to the trigger'])
    wSum(~zSum) = 1;
end
ics = double(1 - cumsum(signMat(:,wIdx)./wSum,2));
% Quartile cuts
% quartCut = exp(-log([4/3, 2, exp(1), 4]))';
% Proportion cuts (5, 25, 50, 63.21, 75, 95)
quartCut = [3.8;3;2;4*exp(-1);1;1/5]/4; 
% Exponential analysis for the auto-correlograms
mdls = zeros(Ns,2); r2 = zeros(Ns,1); qVals = zeros(Ns,length(quartCut));
% Function to fit
sigFun = @(u,v,x) (-1 ./ (1 + exp(-u*(x + v)))) + 1;
ft = fittype(sigFun,'coefficients',{'u','v'},'independent',{'x'});
ftOps = fitoptions(ft); 
ftOps = fitoptions(ftOps,...
    'Lower', [0, -Inf],...
    'Upper', [Inf, 0],...
    'Robust', 'Bisquare',...
    'StartPoint', [1e3, -6e-3],...
    'Display', 'off');
st = fittype('SmoothingSpline'); % Empirical finding of the smoothing parameter
stOps = fitoptions(st); stOps = fitoptions(stOps, 'SmoothingParam',1-10^-9.9);
for cis = 1:Ns
    % Exponential fit for the inverted cumsum
    if zSum(cis)
        
        [stObj, stGof] = fit(shTx, ics(cis,:)', st, stOps);
        r2(cis) = stGof.rsquare;
        smthCumResp = feval(stObj, shTx)';
        quartFlags = smthCumResp >= quartCut;
        quartFlags(~any(quartFlags,2),1) = true;
        %[qSubs, ~] = find(diff(quartFlags'));
        qSubs = arrayfun(@(x) find(diff(quartFlags(x,:)'),1,'first'),...
            (1:size(quartFlags,1))');
        %il = arrayfun(@(x) fit_poly(shTx(x:x+1), ics(cis,x:x+1), 1), qSubs,...
        %    'UniformOutput', 0);il = cat(2,il{:});
        il = arrayfun(@(x) fit_poly(shTx(x:x+1), smthCumResp(x:x+1), 1), qSubs,...
            'UniformOutput', 0);il = cat(2,il{:});
        if ~isempty(il) && size(il,2) == 6 && all(il(1,:) < 0)
            qVals(cis,:) = (quartCut' - il(2,:))./il(1,:);
            continue
        end
        % [fitObj, gof] = fit(shTx, ics(cis,:)', 'exp1', 'Display', 'off');
        [fitObj, gof] = fit(shTx, ics(cis,:)', ft, ftOps);
        mdls(cis,:) = coeffvalues(fitObj); r2(cis) = gof.rsquare;
        % Quartiles cut for exponential distribution (5, 25, 50, 63.21, 75, 95)
        % quartFlags = ics(cis,:) >= quartCut;
        smthCumResp = feval(fitObj, tx)';
        quartFlags = smthCumResp >= quartCut;
        quartFlags(~any(quartFlags,2),1) = true;
        [qSubs, ~] = find(diff(quartFlags'));
        %[~, qSubs] = find(quartFlags,6,'last');
        %il = arrayfun(@(x) fit_poly(shTx(x:x+1), ics(cis,x:x+1), 1), qSubs,...
        %    'UniformOutput', 0);il = cat(2,il{:});
        il = arrayfun(@(x) fit_poly(tx(x:x+1), smthCumResp(x:x+1), 1), qSubs,...
            'UniformOutput', 0);il = cat(2,il{:});
        if isempty(il)
            continue
        end
        qVals(cis,:) = (quartCut' - il(2,:))./il(1,:);
    end
end
qDiff = diff(qVals(:,[1,4]),1,2);
end

% Quartile values extracted from Richard Arnold Johnson; Dean W. Wichern
% (2007). Applied Multivariate Statistical Analysis. Pearson Prentice Hall.
% ISBN 978-0-13-187715-3.