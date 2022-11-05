function [tRes, p] = proportionTest(propCell, p0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ttOpts = {p0, 'Tail', 'right'}; fnOpts = {'UniformOutput', false};
Ncl = size(propCell{1},1);
[tRes, p, ci] = cellfun(@(x) arrayfun(@(y) ttest(double(x(y,:))', ...
    ttOpts{:}), (1:Ncl)',fnOpts{:}), propCell, fnOpts{:});
tRes = repackCells(tRes,5);
p = repackCells(p,5);
    function out = repackCells(cellInput,fl)
        out = cellfun(@(x) cell2mat(x), cellInput, fnOpts{:});
        if fl > 1
            out = cat(2, out{:});
        end
    end
end

