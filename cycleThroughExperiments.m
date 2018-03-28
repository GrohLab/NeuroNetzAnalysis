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
% Eigenvector and Eigenvalue storage:
% 2 dimensions, 6 vectors for three conditions (2 each) and the total
% number of experiments.
vectStack = zeros(2,6,Nex);
Mstack = zeros(7,3,Nex);
MRL = zeros(2,3,Nex);
%% Directories:
rawFiltResponsesPath = fullfile(EphysPath,'continuousEphys');
LFPPath = fullfile(EphysPath,'LFP');
WhiskerPath = fullfile(EphysPath,'Whisker');
AnalysisPath = fullfile(EphysPath,'AnalysisMatFiles');
load([fullfile(EphysPath),'EmptyHead.mat'],'ExampleHead')
for cex = 1:Nex
    fprintf('Dealing with experiment %s...\n',...
        RecDB.Properties.RowNames{cex})
    LFPprobeDepth=ExpDB{{RecDB.AnimalName{cex}}, 'LfpCoord'}(3);
    expDD = discData(cex);
    Nl = round(expDD.LengthInd*(fsLFP/fs));
    if LFPprobeDepth
        ExpName = RecDB.Properties.RowNames{cex};
        %% Loading the raw response and the filtered response
        try
            load([fullfile(rawFiltResponsesPath,ExpName),'.mat'],'fR')
        catch RError
            disp(['No raw response in ',ExpName])
            fR = [];
        end
        try
            load([fullfile(rawFiltResponsesPath,ExpName),'.mat'],'rR')
        catch FError
            disp(['No filtered response in ',ExpName])
            rR = [];
        end
        %% Load LFP and whisker
        [LFP, whisker] = loadLFPAndWhisker(...
            LFPprobeDepth,ExpName,EphysPath);
        
        %% Reconstruct triggers
        lightSignal = getStimPeriods(expDD,fs,fs,'l');
        puffSignal = getStimPeriods(expDD,fs,fs,'p');
        touchSignal = getStimPeriods(expDD,fs,fs,'t');
        whiskingSignal = getStimPeriods(expDD,fs,fs,'w');
        poleSignal = getStimPeriods(expDD,fs,fs,'pl');
        groomingSignal = getStimPeriods(expDD,fs,fs,'g');
        excludeSignal = getStimPeriods(expDD,fs,fs,'x');
        %% Spikes information 'spikesFindingData'
        [~,sp,~] = getStimPeriods(expDD,fs,fs,'s');
        minISI = 1000*(min(diff(sp))/fs); % Minimal ISI in ms.
        ppms = fs/1000;
        %% Building the analysis structures
        spikeFindingData = struct('thresh',[],'minISI',minISI,...
            'spikes',sp,'ppms',ppms,'timestamp',[]);
        filteredResponse = createStructure(fR,fs,ExampleHead);
        RawResponse = createStructure(rR,fs,ExampleHead);
        EEG = createStructure(LFP,fsLFP,ExampleHead);
        Triggers = struct('whisker',whiskingSignal,'light',lightSignal,...
            'puff',puffSignal,'touch',touchSignal,...
            'whisking',whisker,'pole',poleSignal,...
            'grooming',groomingSignal,'exclude',excludeSignal);
        switch RecDB.PhysioNucleus(cex)
            case 'POm'
                NeurType(1,cex) = true;
                notes = {'POm'};
            case 'VPM'
                NeurType(2,cex) = true;
                notes = {'VPM'};
            case 'intermediate'
                NeurType(3,cex) = true;
                notes = {'intermediate'};
            otherwise
                NeurType(3,cex) = true;
                notes = {'other'};
        end 
        % Conditions variable cannot be created yet. The conditions for the
        % awaken state are way complexer than the anaesthesized mice.
        Conditions = {}; 
        save([fullfile(EphysPath,'AnalysisMatFiles',ExpName),'analysis.mat'],...
            'notes','RawResponse','Triggers','filteredResponse','EEG',...
            'spikeFindingData','Conditions')
        continue;
        
        %% Whiskers periods loading
        [whiskPeriods,wp] = getStimPeriods(expDD,fs,fsLFP,'w');
        whiskSP = sp(whiskPeriods(sp));
        nonwhiskSP = sp(~whiskPeriods(sp));
        %% Light Periods loading
        lightDurations = [lightDurations,expDD.LightLength/fs]; %#ok<AGROW>
        %% Puff loading
        
        conditionsIdxs = {sp,whiskSP,nonwhiskSP};
        %% Touch loading
        [touchPeriods, ~, rafTouch] = getStimPeriods(expDD,fs,fsLFP,'t');
        %%
        if RecDB.UsableLFP(cex) %&& RecDB.Light(cex)
            [LFP, whisker] = loadLFPAndWhisker(...
                LFPprobeDepth,ExpName,EphysPath);
            if RecDB.Light(cex)
                [psthStack(:,:,cex),expNLst(cex,:)]=getPSTH(spT,...
                    [lp;expDD.LightLength/fs],whiskPeriods,puffPeriods,LFP,...
                    psthPrev,psthPost,fsLFP);
            end
            [corrSignal(cex,:),~,corrInfo(1,cex),corrInfo(2,cex)] =...
                corrSpLFP(nonwhiskSP,LFP,fsLFP,[0.5,100]);
            % Whisker aligned STAs
            [STAb(cex,:),STAt(cex,:)]=getSTA(wp,LFP,0.01,0.5,fsLFP);
        else
            [~, whisker] = loadLFPAndWhisker(...
                LFPprobeDepth,ExpName,EphysPath);
        end
        filtWhisk = brainwaves(whisker,fsLFP,{'alpha',5,50});
        anaWhisk = hilbert(filtWhisk);
        anaWhiskNorm = anaWhisk./abs(anaWhisk);
        for ccon = 0:2
            MRL(:,ccon+1,cex) = [mean(real(anaWhiskNorm(conditionsIdxs{ccon+1})));...
                mean(imag(anaWhiskNorm(conditionsIdxs{ccon+1})))];
            wrappedSignal = [angle(anaWhisk(conditionsIdxs{ccon+1})),...
                angle(anaWhisk(conditionsIdxs{ccon+1}))-2*pi,...
                angle(anaWhisk(conditionsIdxs{ccon+1}))+2*pi];
            if ccon ~=1
                wrappedSignal = downsample(wrappedSignal,3);
            end
            paraMatrix = emforgmm(wrappedSignal,12,1e-6,1);
            wideAnglePDF = gmmpdf(paraMatrix,[-2*pi,2*pi]);
            magn = wideAnglePDF.pdf/max(wideAnglePDF.pdf);
            angl = wideAnglePDF.getDataDomain;
            eigSignal = complex(cos(angl).*magn, sin(angl).*magn);
            [vectStack(:,ccon*2 +1:(ccon+1)*2,cex),Mstack(:,ccon+1,cex)] =...
                eigenAnalysis(eigSignal,[],false);
            figure;ph=polarhistogram(angle(anaWhisk(conditionsIdxs{ccon+1})));
            hold on;polarplot(wideAnglePDF.getDataDomain,wideAnglePDF.pdf*...
                (max(ph.BinCounts)/max(wideAnglePDF.pdf)))
            figure;
            pn = polarhistogram(angle(anaWhisk),ph.BinEdges);
            pn2 = ph.BinCounts./pn.BinCounts;
            figure;
            polarhistogram('BinEdges',ph.BinEdges,'BinCounts',pn2)
        end
        
        %% Type of cell recorded and its spikes
        switch RecDB.PhysioNucleus(cex)
            case 'POm'
                NeurType(1,cex) = true;
            case 'VPM'
                NeurType(2,cex) = true;
            otherwise
                NeurType(3,cex) = true;
        end
        %         cexNT = NeurName{NeurType(:,cex)};
        %         whiskerPhase = figure('Name',cexNT);
        %         ph = polarhistogram(angle(anaWhisk(conditionsIdxs{2})),'DisplayName','\phi|spike histogram');hold on;
        %         polarplot(atan2(vectStack(2,4,cex),vectStack(1,4,cex)),max(ph.BinCounts),'DisplayName','EigenVector','Marker','o')
        %         polarplot(atan2(Mstack(7,2,cex),Mstack(6,2,cex)),max(ph.BinCounts),'DisplayName','MRL_{pdf}','Marker','o')
        %         polarplot(atan2(MRL(2,2),MRL(1,2)),max(ph.BinCounts),'DisplayName','MRL_{data}','Marker','o')
        %         polarplot(angle(eigSigPerm),abs(eigSigPerm)*max(ph.BinCounts),'DisplayName','GMM fit')
        %         title(sprintf('MRL_{pdf}: %1.4f, L_2/L_1: %1.4f, MRL_{data}: %1.3f',hypot(Mstack(7,2,cex),Mstack(6,2,cex)),...
        %             Mstack(1,2,cex),hypot(MRL(1,2),MRL(2,2))));legend('show')
        %         fn = sprintf('PL_%s_%s',ExpName,NeurName{NeurType(:,cex)});
        %         savefig(whiskerPhase,fullfile('.\Database\WhiskerHistogramFigures\',fn))
        %         close(whiskerPhase)
        %% LFP Analyses
        % Inter spike interval 10 ms
        % [~, ~, spT] = getInitialBurstSpike(sp/fsLFP,0.01);
        
        % Overall, whisking, non whisking
        
        % Spike correlation with the filtered LFP
        
        
        expNLst(cex) = sum(0<lp);
        % PSTH -- account for all the spikes (intra-burst spikes)
        
        
        %             for ws = 1:3
        %                 psthStack(cex,:,ws) = getPSTH(spT,...
        %                     lp,whiskPeriods,LFP,psthPrev,psthPost,fsLFP);
        %             end
    end
    
end

LFPana = struct('CorrelationInformation',struct('AmpAndLoc',corrInfo,...
    'NormalizedCorrelationSignal',corrSignal),...
    'NeuronType',NeurType,...
    'STA',struct('Bursts',STAb,'Tonic',STAt),...
    'PSTH',struct(...
    'Overall',psthStack(1,:,:),...
    'Whisking',psthStack(2,:,:),...
    'NonWhisking',psthStack(3,:,:),...
    'NTrials',expNLst),...
    'EigenAnalysis',struct(...
    'Overall',struct('Vectors',vectStack(:,1:2,:),'Measures',squeeze(Mstack(:,1,:)),'DataMRL',squeeze(MRL(:,1,:))),...
    'Whisking',struct('Vectors',vectStack(:,3:4,:),'Measures',squeeze(Mstack(:,2,:)),'DataMRL',squeeze(MRL(:,2,:))),...
    'NonWhisking',struct('Vectors',vectStack(:,5:6,:),'Measures',squeeze(Mstack(:,3,:)),'DataMRL',squeeze(MRL(:,3,:)))));
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

function [stimPeriods,stimStart,stimRF] = getStimPeriods(dd,fs,fs2,stimString)
% 'Down-sampling' the time indeces for the stimulus/triggers
fact = fs2/fs;
stimStart = [];
stimLength = [];
switch stimString
    case 's'
        stimStart = round(dd.Spikes*fact);
        stimLength = zeros(size(stimStart));
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
    case 'pl'
        if ~isempty(dd.Sections.Pole)
            stimStart = round(dd.Sections.Pole(:,1)*fact);
            stimLength = round((dd.Sections.Pole(:,2)-dd.Sections.Pole(:,1))...
                * fact);
        end
    case 'g'
        if ~isempty(dd.Sections.Grooming)
            stimStart = round(dd.Sections.Grooming(:,1)*fact);
            stimLength = round((dd.Sections.Grooming(:,2)-dd.Sections.Grooming(:,1))...
                * fact);
        end
    case 'x'
        if ~isempty(dd.Sections.Exclude)
            stimStart = round(dd.Sections.Exclude(:,1)*fact);
            stimLength = round((dd.Sections.Exclude(:,2)-dd.Sections.Exclude(:,1))...
                * fact);
        end
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
stimRF = [stimStart',stimLength'];
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
try
    LFP = load([EphysPath, 'LFP\', ExpName, '.mat'],LFP_L5);
    LFP = struct2array(LFP);
catch
    LFP = [];
end
try
    whisker = load([EphysPath, 'Whisker\', ExpName, '.mat'],...
        'WhiskerAngle');
    whisker = struct2array(whisker);
catch
    whisker = [];
end
end

function varStruct = createStructure(chan,fs,exHead)
if ~isempty(chan)
    auxHead = getChannelHeader(exHead,fs,chan);
    varStruct = struct('data',chan,'header',auxHead);
else
    varStruct = struct('data',[],'header',exHead);
end
end

function head = getChannelHeader(exHead,fs,chan)
head = exHead;
head.min = min(chan);
head.max = max(chan);
head.SamplingFrequency = fs;
head.npoints = numel(chan);
end