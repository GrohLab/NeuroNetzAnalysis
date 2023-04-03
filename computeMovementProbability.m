function [outputArg1,outputArg2] = computeMovementProbability(...
    behStack, bvWin, fr, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fnOpts = {'UniformOutput', false};

Nbs = numel(behStack); [Nbt, Ntr] = size(behStack{1}); trigSubs = 1:Ntr;
tmdl = fit_poly([1, Nbt], bvWin, 1); behTx = ((1:Nbt)'.^[1,0])*tmdl;

myRng = @(x) range(x, "all");
mvRng = cellfun(myRng, behStack);

pk_loc = cellfun(@(bs, m) getWaveformCriticalPoints(bs, fr), behStack, fnOpts{:});
pk_loc = cellfun(@(b) cellfun(@(t) t+bvWin(1), b, fnOpts{:}), pk_loc, fnOpts{:});

right_side_flag = cellfun(@(pk) cellfun(@(tr) tr > 0, pk(:,1), ...
    fnOpts{:}), pk_loc, fnOpts{:});
frstPk = cellfun(@(b) cellfun(@(t) t(find(t>0.01, 1, "first")), b(:,1), ...
    fnOpts{:}), pk_loc, fnOpts{:});
frstPk_all = cellfun(@(b) cat(1, b{:}), frstPk, fnOpts{:});
[~, trOrd] = cellfun(@(b) sort(b), frstPk_all, fnOpts{:});

pvpt = cellfun(@(bs, pl, m) arrayfun(@(tr) interp1(behTx, bs(:,tr), ...
    pl{tr,1}, "cubic"), trigSubs, fnOpts{:}), behStack, pk_loc, fnOpts{:});

pdFlag = false(Ntr, Nbs); pmFlag = pdFlag; mvFlag = pdFlag;
pk_dist = cell(Ntr, Nbs); rosFlag = pk_dist; froFlag = pdFlag;
rgFlag = true(size(froFlag));

for cbs = 1:Nbs
    for ctr = trigSubs(:)'
        if ~isempty(right_side_flag{cbs}{ctr}) && ...
                sum(right_side_flag{cbs}{ctr}) > 1
            % sLoc = ctr+zeros(sum(~right_side_flag{cbs}{ctr}),1)-box_sep;
            sVal = pvpt{cbs}{ctr}(~right_side_flag{cbs}{ctr});
            % rLoc = ctr+zeros(sum(right_side_flag{cbs}{ctr}),1)+box_sep;
            rVal = pvpt{cbs}{ctr}(right_side_flag{cbs}{ctr});
            rTime = pk_loc{cbs}{ctr,1}(right_side_flag{cbs}{ctr});
            if ~isempty(sVal) && ~isempty(rVal)
                [~, pmFlag(ctr, cbs)] = ranksum(sVal, rVal);
                sIqr = iqr(sVal); sQs = quantile(sVal,[1,3]/4);
                rosFlag{ctr, cbs} = isWhiskOutlier(rVal, sQs(:), sIqr);
                frPkSubs = find(rTime>0.01,3,"first");
                froFlag(ctr, cbs) = any(rosFlag{ctr,cbs}(frPkSubs));
                pk_dist{ctr, cbs} = distmatrix(sVal, rVal);
                if numel(sVal) > 1 && numel(rVal) > 1
                    pdFlag(ctr, cbs) = ansaribradley(sVal-median(sVal), ...
                        rVal-median(rVal));
                    if froFlag(ctr, cbs)
                        % Is outlier at least a 10th of the signal range?
                        rgFlag(ctr, cbs) = mean(abs((cumsum(any((...
                            pk_dist{ctr, cbs}(:, frPkSubs) > (mvRng(cbs)/10)), ...
                            2)).^[1,0]) * fit_poly([0, size(sVal,1)], [0,1], 1))) ...
                            > 0.5;
                    end
                end
            end
        end
    end
    mvFlag(:, cbs) = (pdFlag(:,cbs) | pmFlag(:,cbs) | froFlag(:,cbs)) & ...
        rgFlag(:, cbs);
end


outputArg1 = inputArg1;
outputArg2 = inputArg2;
end