function [LFPana] = cycleThroughExperiments(ExpDB,RecDB,...
    discData,EphysPath)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
fs = 2e4;
fsLFP = 1e3;

% Number of experiments
[Nex,~] = size(RecDB);
% Recorded neuron for the current experiment (cex). The first row
% corresponds to pom cells, the second to vpm and the third to other.
NeurType = false(3,Nex);NeurName = {'POm','VPM','other'};
corrInfo = zeros(6,Nex);
staWindowSamples = round(0.5*fsLFP);
STAb = zeros(Nex,staWindowSamples);
STAt = STAb;
corrSignal = zeros(Nex,1001);
% Time before the stimulus onset
psthPrev = 0.5;
% Time after the stimulus onset.
psthPost = 1;
% psthStack = zeros(Nex,(psthPrev+psthPost)*fsLFP + 1);
psthStack = zeros(3,(psthPrev+psthPost)*fsLFP + 1,Nex);
expNLst = zeros(Nex,3);
lightDurations = [];
for cex = 1:Nex
    fprintf('Dealing with experiment %s...\n',...
        RecDB.Properties.RowNames{cex})
    if LFPprobeDepth && ~isempty(LFPprobeDepth)
        expDD = discData(cex);
        ExpName = RecDB.Properties.RowNames{cex};
        LFPprobeDepth=ExpDB{{RecDB.AnimalName{cex}}, 'LfpCoord'}(3);
        sp = round(expDD.Spikes*(fsLFP/fs));
        spT = false(1,Nl);
        spT(sp) = true;
        %% Whiskers periods loading
        [whiskPeriods,wp] = getStimPeriods(expDD,fs,fsLFP,'w');
        whiskSP = sp(whiskPeriods(sp));
        nonwhiskSP = sp(~whiskPeriods(sp));
        %% Light Periods loading
        [lightPeriods,lp] = getStimPeriods(expDD,fs,fsLFP,'l'); %#ok<ASGLU>
        lightDurations = [lightDurations,expDD.LightLength/fs]; %#ok<AGROW>
        %% Puff loading
        [puffPeriods, ~] = getStimPeriods(expDD,fs,fsLFP,'p');
        if RecDB.UsableLFP(cex) %&& RecDB.Light(cex)
            [LFP, whisker] = loadLFPAndWhisker(...
                LFPprobeDepth,ExpName,EphysPath);
            Nl = length(LFP);
            [psthStack(:,:,cex),expNLst(cex,:)]=getPSTH(spT,...
                [lp;expDD.LightLength/fs],whiskPeriods,puffPeriods,LFP,...
                psthPrev,psthPost,fsLFP);
            [corrSignal(cex,:),~,corrInfo(1,cex),corrInfo(2,cex)] =...
                corrSpLFP(nonwhiskSP,LFP,fsLFP,[0.5,100]);
            % Whisker aligned STAs
            [STAb(cex,:),STAt(cex,:)]=getSTA(wp,LFP,0.01,0.5,fsLFP);
        else
            [~, whisker] = loadLFPAndWhisker(...
                LFPprobeDepth,ExpName,EphysPath);
            Nl = length(whisker);
        end
        []
        %% Type of cell recorded and its spikes
        switch RecDB.PhysioNucleus(cex)
            case 'POm'
                NeurType(1,cex) = true;
            case 'VPM'
                NeurType(2,cex) = true;
            otherwise
                NeurType(3,cex) = true;
        end
        cexNT = NeurName{NeurType(:,cex)};
        % Translating the spike times to 1 kHz (or fsLFP)
        
        %% LFP Analyses
        % Inter spike interval 10 ms
        % [~, ~, spT] = getInitialBurstSpike(sp/fsLFP,0.01);
        
        % Overall, whisking, non whisking
        conditionsIdxs = {sp,whiskSP,nonwhiskSP};
        % Spike correlation with the filtered LFP
        
        
        expNLst(cex) = sum(0<lp);
        % PSTH -- account for all the spikes (intra-burst spikes)
        
        
        %             for ws = 1:3
        %                 psthStack(cex,:,ws) = getPSTH(spT,...
        %                     lp,whiskPeriods,LFP,psthPrev,psthPost,fsLFP);
        %             end
    end
    
end
% LFPana = struct('CorrelationInformation',struct('AmpAndLoc',corrInfo,...
%     'NormalizedCorrelationSignal',corrSignal),...
%     'NeuronType',NeurType,...
%     'STA',struct('Bursts',STAb,'Tonic',STAt),...
%     'PSTH',struct('Overall',squeeze(psthStack(:,:,1)),...
%                   'Whisking',squeeze(psthStack(:,:,2)),...
%                   'NonWhisking',squeeze(psthStack(:,:,3)),...
%                   'NTrials',expNLst));
LFPana = struct('CorrelationInformation',struct('AmpAndLoc',corrInfo,...
    'NormalizedCorrelationSignal',corrSignal),...
    'NeuronType',NeurType,...
    'STA',struct('Bursts',STAb,'Tonic',STAt),...
    'PSTH',struct('Overall',psthStack(1,:,:),...
    'Whisking',psthStack(2,:,:),...
    'NonWhisking',psthStack(3,:,:),...
    'NTrials',expNLst));
% figure;histogram(lightDurations);
end

function [fk, k, amp, pos] = corrSpLFP(sp,LFP,fsLFP,freqBand)
LFPf = brainwaves(LFP,fsLFP,{'alpha',freqBand(1),freqBand(2)});
LFPf = zscore(LFPf);
spT = zeros(1,length(LFP));
spT(sp) = true;
[fk, k] = xcorr(spT,LFPf,'coef');
% Dirty approach --> Cut off the signal
fk = fk(k >= -500 & k <= 500);
k = k(k >= -500 & k <= 500);
[amp,lg] = max(fk);
pos = k(lg);
end

function [stimPeriods,stimStart] = getStimPeriods(dd,fs,fs2,stimString)
fact = fs2/fs;
switch stimString
    case 'w'
        stimStart = round(dd.WhiskingStart*fact);
        stimLength = round(dd.WhiskingLength*fact);
    case 'l'
        stimStart = round(dd.LightStart*fact);
        stimLength = round(dd.LightLength*fact);
    case 'p'
        stimStart = round(dd.PuffStart*fact);
        stimLength = round(dd.PuffLength*fact);
    case 't'
        stimStart = round(dd.TouchStart*fact);
        stimLength = round(dd.TouchLength*fact);
    otherwise
        stimPeriods = dd.LengthTime;
        stimStart = 0;
        fprintf('No stimuli recognized...')
        return;
end
stimPeriods = false(1,round(dd.LengthTime*fs2));
for counter = 1:length(stimStart)
    stimPeriods(...
        stimStart(counter):stimStart(counter)+stimLength(counter)) = true;
end
end


function [LFP,whisker] = loadLFPAndWhisker(LFPprobeDepth,ExpName,EphysPath)
% LFP is recorded in 16 linear channels along cortex. The order from white
% matter to pia is:
LFPsort=[6, 11, 3, 14, 1, 16, 2, 15, 5, 12, 4, 13, 7, 10, 8, 9];
% The individual channels depth is:
LFPdepth=LFPprobeDepth - (50:100:1550);
% Consequently for ~L5 lfp at 850 ?m depth, when LFPprobeDepth
% is 1600 would be:
LFPchIdx = LFPsort(find(LFPdepth < 900, 1));
% Loading only the L5 LFP from the 16 channels file
LFP_L5 = sprintf('LFPch%d',LFPchIdx);
LFP = load([EphysPath, 'LFP\', ExpName, '.mat'],LFP_L5);
whisker = load([EphysPath, 'Whisker\', ExpName, '.mat'],...
    'WhiskerAngle');
whisker = struct2array(whisker);
LFP = struct2array(LFP);
end