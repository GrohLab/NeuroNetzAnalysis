function [H, P] = statTests(dStack, condFlag, timeFlags, varargin)
%STATTESTS function description
% Emilio Isaias-Camacho @ GrohLab
H = [];
P = [];
%% Validating the inputs
p = inputParser;

defaultTest = 'kstest';
validTests = {'kstest','mcnemar','chi2','binomial','wilcoxon'};
checkTest = @(x) any(validatestring(x,validTests));

checkStack = @(x) any([islogical(x), size(x) >= 1, numel(size(x)) == 3]);

[~, Nt, NTa] = size(dStack);

checkCondFlag = @(x) all([size(x,1) == NTa, islogical(x)]);
checkTimeFlags = @(x) all([size(x,1) == 2,size(x,2) == Nt, islogical(x)]);


addRequired(p, 'dStack', checkStack)
addRequired(p, 'condFlag', checkCondFlag)
addRequired(p, 'timeFlags', checkTimeFlags)
addOptional(p, 'test', defaultTest, checkTest)

p.KeepUnmatched = true;

parse(p,dStack, condFlag, timeFlags, varargin{:})
testType = p.Results.test;

%% Preparatory variables


% Number of alignment points per conditions
Na = sum(condFlag, 1);

% Test selector
switch testType
    case validTests{1} % KStest
        statFun = @kstest2;
    case validTests{2} % McNemar
        fprintf(1,'Not yet implemented. Selecting kstest\n')
        statFun = @kstest2;
        
    case validTests{3} % Chi^2
        
    case validTests{4} % Binomial
        
    case validTests{5} % Wilcoxon
        
end

% Output variable

% Number of conditions to cycle through
Ncond = size(condFlag, 2);
for cc1 = 1:Ncond
    sponCountA =...
        squeeze(sum(dStack(2:end,timeFlags(1,:),condFlag(:,cc1)),2));
    evokCountA =...
        squeeze(sum(dStack(2:end,timeFlags(2,:),condFlag(:,cc1)),2));
    cc2 = cc1 + 1;
    while cc2 <= Ncond
        % Comparing condition A versus condition B; spontaneous and evoked
        sponCountB =...
            squeeze(sum(dStack(2:end,timeFlags(1,:),condFlag(:,cc2)),2));
        evokCountB =...
            squeeze(sum(dStack(2:end,timeFlags(2,:),condFlag(:,cc2)),2));
        
    end
    % Same condition testing
    [P, H] = signrank(sponCountA, evokCountA);
    shufIdx = randperm()
    [P, H] = signrank(SponCountA, evokCountA);
end
end

