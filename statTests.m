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
checkTimeFlags = @(x) all([size(x,1) == 2,size(x,2) == Nt, islogical(x)]);


addRequired(p, 'dStack', checkStack)
addRequired(p, 'condFlag', checkCondFlag)
addRequired(p, 'timeFlags', checkTimeFlags)
addOptional(p, 'test', defaultTest, checkTest)

p.KeepUnmatched = true;

parse(p,dStack, condFlag, timeFlags, varargin{:})
testType = p.Results.test;

%% Preparatory variables
% Number of conditions to cycle through
Ncond = size(condFlag, 2);
for cc1 = 1:Ncond
    for cc2 = cc1:Ncond
        if cc1 == cc2
            % Comparing spontaneous versus evoked from the same condition
            sponCount =...
                squeeze(sum(dStack(2:end,timeFlags(1,:),condFlag(:,cc1)),2));
            evokCount =...
                squeeze(sum(dStack(2:end,timeFlags(2,:),condFlag(:,cc1)),2));
            
            continue
        end
        % Comparing condition A versus condition B; spontaneous and evoked
    end
end
end

