function [PSTHStruct] = getPopPSTH(dSt, condSt, configStruct)
%GETPOPPSTH returns the peri-stimulus triggered histogram for all the
%experiments collected in the population stack and conditions it according
%to both the configuration and conditioning structure.
%   Emilio Isaías-Camacho @GrohLab 2019

%% Initialization
% Gathering the necessary information to compute the PSTH combinations
% according to the user selection.
[Nv, Nts, Nap] = size(dSt.Stack);
fs = configStruct.SamplingFrequencies(1); 
Ncv = size(condSt.CVDFlags,1); % Number of conditioning variables
finishedButton = false;
% To have a pure triggered PSTH (only taking random occurances of the
% ignored variables) we need to take those trials which were NOT excluded
% AND which had NO conditioning variables happening during the selected
% window of interest. De Morgan law: ~a & ~b = ~(a | b)
cleanIdx = ~(condSt.ExcludeFlags | any(condSt.CVDFlags,1));
cl_swaps = sum(cleanIdx);
fprintf('Percentage of ''clean'' trials: %.2f%% (%d/%d)\n',100*cl_swaps/Nap,...
    cl_swaps, Nap)
% cleanPSTH = sum(dSt.Stack(:,:,cleanIdx),3)./cl_swaps;
tx = (0:Nts-1) * 1/fs - configStruct.ViewWindow(1);
%% Possible combinations
% $sum_{i=0}^{N_{cv}-1}(2^i)+1$
possCombi = sum(2.^(0:Ncv-1)) + 1;
%% Computing the conditional PSTHs with user dependent selection.
confCell = squeeze(struct2cell(configStruct.ConditionWindow));
Nexp = numel(condSt.ExperimentCuts) - 1;
currentPSTH = zeros(Nv, Nts, Nexp, 'single');
auxResPSTH = struct('ConditionCombination',{{'Control'}},...
    'PSTHstack',currentPSTH,...
    'TrialsPerExperiment', condSt.ExperimentCuts(2:Nexp+1));
auxResPSTH = repmat(auxResPSTH,1,possCombi);
auxCnt = 1;
% Clean PSTH stack
for cexp = 1:Nexp
        curExp = sum(condSt.ExperimentCuts(1:cexp))+1:...
            sum(condSt.ExperimentCuts(1:cexp+1));
        currentPSTH(:,:,cexp) = sum(dSt.Stack(:,:,curExp(cleanIdx(curExp))),3);
        auxResPSTH(auxCnt).TrialsPerExperiment(cexp) =...
            sum(cleanIdx(curExp));
end
auxResPSTH(auxCnt).PSTHstack = currentPSTH;
auxCnt = auxCnt + 1;
while ~finishedButton && auxCnt <= possCombi
    % While the possible combinations are not exhausted, the loop continues
    % to prompt the user for interesting combinations. A future plan is to
    % compute all combinations and then let the user choose which
    % combinations to plot. That seems like a good plan.
    scvFlags = true(1,Ncv);
    if Ncv > 1
        [combSubs, iOk] = listdlg('ListString',confCell(1,:)',...
            'PromptString','Select a conditioning variable combination',...
            'SelectionMode','multiple',...
            'CancelString','Finish',...
            'OKString','Create',...
            'Name','Combination',...
            'ListSize',[160,15*Ncv],...
            'InitialValue',1);
        scvFlags(combSubs) = false;
        if iOk && sum(~scvFlags)
            scvFlags = ~scvFlags;
        else
            finishedButton = true;
            continue
        end
    end
    inTrials = ~condSt.ExcludeFlags & any(condSt.CVDFlags(scvFlags,:),1);
    % Single experiment PSTH
    currentPSTH = zeros(Nv, Nts, Nexp, 'single');
    for cexp = 1:Nexp
        curExp = sum(condSt.ExperimentCuts(1:cexp))+1:...
            sum(condSt.ExperimentCuts(1:cexp+1));
        currentPSTH(:,:,cexp) = sum(dSt.Stack(:,:,curExp(inTrials(curExp))),3);%/...
            % sum(inTrials(curExp));
        auxResPSTH(auxCnt).TrialsPerExperiment(cexp) =...
            sum(inTrials(curExp));
    end
    % currentPSTH = sum(dSt.Stack(:,:,inTrials),3)/sum(inTrials);
    % The population stack can be performed easily outside the function.
    auxResPSTH(auxCnt).PSTHstack = currentPSTH;
    auxResPSTH(auxCnt).ConditionCombination =...
        {confCell{[true,false],scvFlags}};
    auxCnt = auxCnt + 1;
end
%% Wrapping up the results
% The empty allocated spaces will be deleted.
cres = possCombi;
emptyFlag = true;
while cres > 1 && emptyFlag
    emptyFlag = sum(strcmpi(auxResPSTH(cres).ConditionCombination,...
        'none'));
    try
        auxResPSTH(cres*emptyFlag) = [];
        cres = cres - 1;
    catch
        fprintf('%d combinations unsued\n',possCombi-cres)
    end
end
PSTHStruct = struct('BinSize',configStruct.BinSize,...
    'TimeAxis',tx,...
    'Trigger',configStruct.Trigger.Name,...
    'PSTHs',auxResPSTH,...
    'SignalIDs',{dSt.SignalIDs});
end