function [PSTHStruct] = getPopPSTH(dSt,condSt,configStruct)
%GETPOPPSTH returns the peri-stimulus triggered histogram for all the
%experiments collected in the population stack and conditions it according
%to both the configuration and conditioning structure.
%   Emilio Isaías-Camacho @GrohLab 2019

%% Initialization
% Gathering the necessary information to compute the PSTH combinations
% according to the user selection.
[~, Nts, Nap] = size(dSt.Stack);
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
cleanPSTH = sum(dSt.Stack(:,:,cleanIdx),3)./cl_swaps;
tx = (0:Nts-1) * 1/fs - configStruct.ViewWindow(1);
possCombi = sum(2.^(0:Ncv-1)) + 1;
auxResPSTH = struct('ConditionCombination','None','PSTH',cleanPSTH);
auxResPSTH = repmat(auxResPSTH,1,possCombi);
auxCnt = 2;
confCell = squeeze(struct2cell(configStruct.ConditionWindow));
%% Computing the conditional PSTHs with user dependent selection.
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
    currentPSTH = sum(dSt.Stack(:,:,inTrials),3)/sum(inTrials);
    auxResPSTH(auxCnt).PSTH = currentPSTH;
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
    'PSTHs',auxResPSTH);
end

% binWidth = ceil(configStruct.BinSize * fs);
% Nbin = ceil(Nts/binWidth);