function [condFlags, Na] = conditionBlock(condStruct, idx)
%CONDITIONBLOCK looks for time gaps in which other stimulus falls into. For
%example, separating whisker-control, laser induction, and whisker-read-out
%with a logical array.
%   [condFlags, Na] = conditionBlock(condStruct, idx)
NTa = size(condStruct(idx).Triggers,1);
[biGaps, gapSubs] = sort(diff(condStruct(idx).Triggers(:,1)), 'descend');
gapFlag = abs(zscore(biGaps)) > 3;
mainGaps = sort(gapSubs(gapFlag), 'ascend');
Ng = numel(mainGaps);
Ncond = Ng+1;
condFlags = false(NTa, Ncond);
initSub = 1;
for cg = 1:Ng
    condFlags(initSub:mainGaps(cg),cg) = true;
    initSub = mainGaps(cg) + 1;
end
condFlags(initSub:NTa, Ncond) = true;
Na = sum(condFlags,1);
end

