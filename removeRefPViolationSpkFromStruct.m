function [cleanSpkStr] = removeRefPViolationSpkFromStruct(spkStruct, varargin)
%REMOVEREFVIOLATIONSPKFROMSTRUCT cleans all the spike trains in the
%relative spike structure which are below $rp$. 
%   cleanSpkStr = removeRefViolationSpkFromStruct(spkStruct);
%   cleanSpkStr = removeRefViolationSpkFromStruct(spkStruct, 'RefPeriod', rp);
%Emilio Isaias-Camacho@GrohLab 2022
%% Input validation
p = inputParser;
m =  1e-3;

fldNames = {'name', 'SpikeTimes'};
fnOpts = {'UniformOutput', false};

defRP = m;
checkRP = @(x) all([isnumeric(x), x>0, numel(x) == 1]);
checkSpkStruct = @(x) all([arrayfun(@(y) isstruct(y), x), ...
    arrayfun(@(y) all(isfield(y, fldNames)), x)]);

addRequired(p, 'spkStruct', checkSpkStruct)
addParameter(p, 'RefPeriod', defRP, checkRP)

parse(p, spkStruct, varargin{:})

relSpkTmsStruct = p.Results.spkStruct;
rp = p.Results.RefPeriod;
%% Spike cleaning
% Inter-spike intervals (ISIs)
isi = arrayfun(@(x) cellfun(@(y) diff(y), x.SpikeTimes, fnOpts{:}), ...
    relSpkTmsStruct, fnOpts{:});
% Comparing the ISIs with the chosen refractory period.
rpvFlag = cellfun(@(x) cellfun(@(y) y < rp, x, fnOpts{:}), isi, fnOpts{:});
% Removing the spikes
cleanSpkTms = cellfun(@(x,y) cellfun(@(u,v) u([true&numel(u), ~v]), ...
    x, y, fnOpts{:}), {relSpkTmsStruct.SpikeTimes}, rpvFlag, ...
    fnOpts{:}); cleanSpkStr = relSpkTmsStruct;
% Organising the output
for ci = 1:numel(relSpkTmsStruct)
    cleanSpkStr(ci).SpikeTimes = cleanSpkTms{ci};
end

end