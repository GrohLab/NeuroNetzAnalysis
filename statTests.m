function [H, P] = statTests(dStack, condFlag, timeFlags, varargin)
%STATTESTS function description
% Emilio Isaias-Camacho @ GrohLab
H = [];
P = [];
%% Validating the inputs
p = inputParser;

defaultTest = 'kstest';
validTests = {'kstest','mcnemar','chi2','binomial'};
checkTest = @(x) any(validatestring(x,validTests));

checkStack = @(x) any([islogical(x), size(x) >= 1, numel(size(x)) == 3]);

[~, Nt, NTa] = size(dStack);

checkCondFlag = @(x) all([size(x,1) == NTa, islogical(x)]);
checkTimeFlags = @(x) all([length(x) == Nt, islogical(x)]);


addRequired(p, 'dStack', checkStack)
addRequired(p, 'condFlag', checkCondFlag)
addRequired(p, 'timeFlags', checkTimeFlags)
addOptional(p, 'test', defaultTest, checkTest)

p.KeepUnmatched = true;

parse(p,dStack, condFlag, timeFlags, varargin{:})
testType = p.Results.test;

%% Preparatory variables

end

