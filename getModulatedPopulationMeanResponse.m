function popResponse = getModulatedPopulationMeanResponse(spkCount, rmdFlag, varargin)
%getModulatedPopulationMeanResponse uses the counts and flags to retrieve
%the modulated clusters from the given counts.
%       popResponse = getModulatedPopulationMeanResponse(spkCount, rmdFlag)
%   INPUTS:
%       - spkCount - CxE cell matrix the spike counts within E certain time
%                    windows with C conditions. 
%       - rmdFlag - Nx3 logic array specifying the neurons which are
%                   responsive (:,1), significantly modulated (:,2), and
%                   either positively or negatively modulated.
%    OUTPUTS:
%       - popResponse - structure containing the mean and confidence limits
%                       for the given spike counts.
% Emilio Isaias-Camacho @GrohLab 2021

%% Parse inputs
p = inputParser;

checkFlags = @(x) all([isnumeric(x) | islogical(x), size(x,2) == 3]);


p.addRequired('spkCount',@iscell);
p.addRequired('rmdFlag',checkFlags);

p.parse(spkCount, rmdFlag);

spkCount = p.Results.spkCount;
rmdFlag = p.Results.rmdFlag;
%% Validation for agrument interaction
fnOpts = {'UniformOutput',false}; 
popResponse = struct('Mean',[],'Confidence',[]);
countMatSzs = cellfun(@(x) size(x,1), spkCount);
if any((countMatSzs - countMatSzs) ~= 0)
    fprintf(2,'The given cluster counts contain different cluster number!\n')    
    qstAns = questdlg('Continue?',[],'Yes','No','No');
    if strcmp(qstAns,'No') 
        return;
    end
end

if any((countMatSzs - size(rmdFlag,1)) ~= 0)
    fprintf(2,'Impossible to continue!\n')
    fprintf(1,'The given flags contain logical values outside some count')
    fprintf(1,' matrices.\n')
    return;
end
%% Stimulus response evolution
modSpks = arrayfun(@(x)...
    cellfun(@(y) y(all(xor(rmdFlag,[0,0,x]),2),:),...
    spkCount(:,2), fnOpts{:}), [false,true], fnOpts{:}); 
modSpks = cat(2, modSpks{:}); Ntrials = cellfun(@(x) size(x,2), modSpks, fnOpts{:});
poissonDistributions =...
    cellfun(@(x,u)...
    arrayfun(@(y) fitdist(x(:,y), 'Poisson'), 1:u, fnOpts{:}),...
    modSpks, Ntrials, fnOpts{:});
poissonDistributions = cellfun(@(x) cat(2, x{:}), poissonDistributions, fnOpts{:});
rlambda = cellfun(@(x) cat(2, x.lambda), poissonDistributions, fnOpts{:});
rconfidence = cellfun(@(x)...
    arrayfun(@(y) y.paramci, x, fnOpts{:}),...
    poissonDistributions, fnOpts{:}); 
rconfidence = cellfun(@(x) cat(2, x{:}), rconfidence, fnOpts{:});
ffact = cellfun(@(x) var(x,'omitnan')./mean(x,'omitnan'), spkCount, fnOpts{:});
% Standard error of the mean
sem = cellfun(@(x,y) sqrt(x./size(y,1)), rlambda, modSpks, fnOpts{:});
popResponse = struct('Mean', rlambda, 'Confidence', rconfidence,...
    'FanoFactor',ffact,'SEM', sem);
end

