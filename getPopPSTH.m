function [psthStr] = getPopPSTH(dSt,condSt,configStruct)
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
cleanPSTH = sum(dSt.Stack(:,:,cleanIdx),3)./cl_swaps;
tx = (0:Nts-1) * 1/fs - configStruct.ViewWindow(1);
binWidth = ceil(configStruct.BinSize * fs);
Nbin = ceil(Nts/binWidth);
btx = (0:Nbin-1) * configStruct.BinSize - configStruct.ViewWindow(1);
possCombi = sum(2.^(0:Ncv-1)) + 1;
binCleanPSTH = binTime(cleanPSTH,configStruct.BinSize,fs);
auxResPSTH = struct('ConditionCombination','None','PSTH',binCleanPSTH);
auxResPSTH = repmat(auxResPSTH,1,possCombi);
auxCnt = 2;
confCell = struct2cell(configStruct.ConditionWindow);
while ~finishedButton && auxCnt <= possCombi
    scvFlags = true(1,Ncv);
    if Ncv > 1
        [combSubs, iOk] = listdlg('ListString',confCell{1,:},...
            'PromtString','Select a conditioning variable combination',...
            'SelectionMode','multiple',...
            'CancelString','Finish',...
            'OKString','Create',...
            'Name','Combination');
        scvFlags(combSubs) = false;
        if iOk && sum(~scvFlags)
            scvFlags = ~scvFlags;
        else
            finishedButton = true;
            continue
        end
    end
    inTrials = ~condSt.ExcludeFlags & any(condSt.CVDFlags(scvFlags,:),1);
    uniquePSTH = sum(dSt.Stack(:,:,inTrials),3)/sum(inTrials);
    auxResPSTH(auxCnt).PSTH =...
        binTime(uniquePSTH,configStruct.BinSize,fs);
    auxResPSTH(auxCnt).ConditionCombination =...
        confCell{[true,false],scvFlags};
    auxCnt = auxCnt + 1;
end

psthStr = struct('BinSize',configStruct.BinSize,...
    'TimeAxes',struct('Trigger',tx,'PSTH',btx),...
    'Trigger',cleanPSTH(1,:));
end

