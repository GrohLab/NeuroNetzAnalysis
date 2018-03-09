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
psthStack = zeros(2,(psthPrev+psthPost)*fsLFP + 1,Nex);
expNLst = zeros(Nex,2);
for cex = 1:Nex
    fprintf('Dealing with experiment %s...\n',...
        RecDB.Properties.RowNames{cex})
    if RecDB.UsableLFP(cex) && RecDB.Light(cex)
        LFPprobeDepth=ExpDB{{RecDB.AnimalName{cex}}, 'LfpCoord'}(3);
        if LFPprobeDepth && ~isempty(LFPprobeDepth)
            expDD = discData(cex);
            ExpName = RecDB.Properties.RowNames{cex};
            [LFP, whisker] = loadLFPAndWhisker(...
                LFPprobeDepth,ExpName,EphysPath);
            Nl = length(LFP);
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
            sp = round(expDD.Spikes*(fsLFP/fs));
            %% Whiskers periods loading
            [whiskPeriods,wp] = getWhiskPeriods(expDD,fs,fsLFP,Nl);
            whiskSP = sp(whiskPeriods(sp));
            nonwhiskSP = sp(~whiskPeriods(sp));
            %% LFP Analyses
            % Inter spike interval 10 ms
            [~, ~, spT] = getInitialBurstSpike(sp/fsLFP,0.01);
            % Overall, whisking, non whisking
            conditionsIdxs = {sp,whiskSP,nonwhiskSP};
            % Spike correlation with the filtered LFP
            [corrSignal(cex,:),~,corrInfo(1,cex),corrInfo(2,cex)] =...
                corrSpLFP(nonwhiskSP,LFP,fsLFP,[0.5,100]);
            % Whisker aligned STAs
            [STAb(cex,:),STAt(cex,:)]=getSTA(wp,LFP,0.01,0.5,fsLFP);
            [~,lp] = getLightPeriods(expDD,fs,fsLFP,Nl);
            expNLst(cex) = sum(0<lp);
            [psthStack(:,:,cex),expNLst(cex,1),expNLst(cex,2)]=getPSTH(spT,...
                lp,whiskPeriods,LFP,psthPrev,psthPost,fsLFP);
%             for ws = 1:3
%                 psthStack(cex,:,ws) = getPSTH(spT,...
%                     lp,whiskPeriods,LFP,psthPrev,psthPost,fsLFP);
%             end
        end
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
    'NTrials',expNLst));
%                   'Whisking',psthStack(1,:,:),...
%                   'NonWhisking',psthStack(2,:,:),...
                  

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

function [lightPeriods,lp] = getLightPeriods(dd,fs,fs2,N)
lp = round(dd.LightStart*(fs2/fs));
wl = round(dd.LightLength*(fs2/fs));
lightPeriods = false(1,N);
for counter = 1:length(lp)
    lightPeriods(lp(counter):lp(counter)+wl(counter)) = true;
end
end

function [whiskPeriods,wp] = getWhiskPeriods(dd,fs,fs2,N)
wp = round(dd.WhiskingStart*(fs2/fs));
wl = round(dd.WhiskingLength*(fs2/fs));
whiskPeriods = false(1,N);
for counter = 1:length(wp)
    whiskPeriods(wp(counter):wp(counter)+wl(counter)) =...
        true(1,wl(counter)+1);
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