function [trialSubs, chTrilSubs] = equiliseTrials(trialFlags, excludeFlags)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fnOpts = {'UniformOutput', false};

if ~exist('excludeFlags', 'var') || isempty(excludeFlags) ||...
        size(excludeFlags, 2) ~= 1
    excludeFlags = false(size(trialFlags,1),1);
end
checkTF = @(x) all([islogical(x), ismatrix(x) | isvector(x)]);
if ~checkTF(trialFlags)
    error("Trial flags must be a NxC logical matrix! N - Trials, C - conditions")    
end
[Nt, Nc] = size(trialFlags); Nte = size(excludeFlags, 1);
if Nt ~= Nte
    error("Trial flags and exclude flags must have the same row number!")
end

effFlags = arrayfun(@(x) ~excludeFlags & trialFlags(:,x), 1:Nc, fnOpts{:});
effFlags = cat(2, effFlags{:}); Na = sum(effFlags,1); Nma = min(Na);

trialID = arrayfun(@(x) find(~excludeFlags & trialFlags(:,x)), 1:Nc, ...
    fnOpts{:});
chTrilSubs = arrayfun(@(x) sort(randperm(x, Nma), 'ascend')', Na, ...
    fnOpts{:}); 
trialSubs = cellfun(@(x,y) x(y), trialID, chTrilSubs, fnOpts{:}); 
trialSubs = cat(2, trialSubs{:});

end