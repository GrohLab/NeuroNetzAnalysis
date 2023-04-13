 function [movStruct] = computeMovementProbability(behStack, pairStim, ...
     bvWin, fr, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fnOpts = {'UniformOutput', false};
isWhiskOutlier = @(x, Qs, Qr) (Qs(1,:) - 1.5*Qr) > x(:)' | ...
    (Qs(2,:) + 1.5*Qr) < x(:)';
sqzePge = @(mat, p) squeeze(mat(:,:,p));
pge2cell = @(mat) arrayfun(@(pp) sqzePge(mat, pp), (1:size(mat,3))', fnOpts{:});
[trigSubs, condIDs] = find(pairStim);
condID = unique(condIDs); Nccond = numel(condID);

Nbs = numel(behStack); [Nbt, Ntr] = size(behStack{1}); 
tmdl = fit_poly([1, Nbt], bvWin, 1); behTx = ((1:Nbt)'.^[1,0])*tmdl;

myRng = @(x) range(x, "all");
mvRng = cellfun(myRng, behStack);

pk_loc = cellfun(@(bs, m) getWaveformCriticalPoints(bs, fr), behStack, fnOpts{:});
pk_loc = cellfun(@(b) cellfun(@(t) t+bvWin(1), b, fnOpts{:}), pk_loc, fnOpts{:});

pvpt = cellfun(@(bs, pl, m) arrayfun(@(tr) interp1(behTx, bs(:,tr), ...
    pl{tr,1}, "cubic"), 1:size(pl,1), fnOpts{:}), behStack, pk_loc, fnOpts{:});

right_side_flag = cellfun(@(pk) cellfun(@(tr) tr > 0, pk(:,1), ...
    fnOpts{:}), pk_loc, fnOpts{:});
frstPk = cellfun(@(b) cellfun(@(t) t(find(t > 0.01, 1, "first")), b(:,1), ...
    fnOpts{:}), pk_loc, fnOpts{:});
frstPk_all = cellfun(@(b) cat(1, b{:}), frstPk, fnOpts{:});
[~, trOrd] = cellfun(@(b) sort(b), frstPk_all, fnOpts{:});

pdFlag = false(Ntr, Nccond, Nbs); pmFlag = pdFlag; froFlag = pdFlag;
pk_dist = cell(Ntr, Nccond, Nbs); rosFlag = pk_dist; 
rgFlag = true(size(froFlag));

for cbs = 1:Nbs
    for ccond = condID(:)'
        for ctr = reshape(trigSubs(condIDs==ccond), 1, [])
            if ~isempty(right_side_flag{cbs}{ctr}) && ...
                    sum(right_side_flag{cbs}{ctr}) > 1
                % sLoc = ctr+zeros(sum(~right_side_flag{cbs}{ctr}),1)-box_sep;
                sVal = pvpt{cbs}{ctr}(~right_side_flag{cbs}{ctr});
                % rLoc = ctr+zeros(sum(right_side_flag{cbs}{ctr}),1)+box_sep;
                rVal = pvpt{cbs}{ctr}(right_side_flag{cbs}{ctr});
                rTime = pk_loc{cbs}{ctr,1}(right_side_flag{cbs}{ctr});
                if ~isempty(sVal) && ~isempty(rVal)
                    [~, pmFlag(ctr, ccond, cbs)] = ranksum(sVal, rVal);
                    sIqr = iqr(sVal); sQs = quantile(sVal,[1,3]/4);
                    rosFlag{ctr, ccond, cbs} = isWhiskOutlier(rVal, sQs(:), sIqr);
                    frPkSubs = find(rTime>0.01,3,"first");
                    froFlag(ctr, ccond, cbs) = ...
                        any(rosFlag{ctr, ccond, cbs}(frPkSubs));
                    pk_dist{ctr, ccond, cbs} = distmatrix(sVal, rVal);
                    if numel(sVal) > 1 && numel(rVal) > 1
                        pdFlag(ctr, ccond, cbs) = ansaribradley(sVal-median(sVal), ...
                            rVal-median(rVal));
                        if froFlag(ctr, ccond, cbs)
                            % Is outlier at least a 10th of the signal range?
                            rgFlag(ctr, ccond, cbs) = mean(abs((cumsum(any((...
                                pk_dist{ctr, ccond, cbs}(:, frPkSubs) > (mvRng(cbs)/10)), ...
                                2)).^[1,0]) * fit_poly([0, size(sVal,1)], [0,1], 1))) ...
                                > 0.5;
                        end
                    end
                end
            end
        end
    end
end

mvFlag = (pdFlag | pmFlag | froFlag) & rgFlag;
movStruct = struct('MovmentFlags', pge2cell(mvFlag), ...
    'FirstMovement', frstPk_all, 'TrialOrder', trOrd, ...
    'MedianFlag', pge2cell(pmFlag), 'DispersionFlag', pge2cell(pdFlag), ...
    'PeakOutlierFlag', pge2cell(froFlag), 'RangeFlag', pge2cell(rgFlag));
end