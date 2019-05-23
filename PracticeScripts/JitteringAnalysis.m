%% Constant values
% Prefixes
k = 1e3; % Kilo
m = 1e-3; % milli
zs = @(x)mean(x)/std(x);
addFst = @(x,y)cat(find(size(x)~=1),y,x);
addLst = @(x,y)cat(find(size(x)~=1),x,y);
uninterestingVals = @(x)~(isempty(x) || isnan(x) || isinf(x) || x == 0);
%% Select data file
defaultPath = 'E:\Data\Jittering'; % Path where you would like to open a windows explorer window to select your data file
[fileName, filePath] =...
    uigetfile({'*.mat';'*.smr';'*.smrx'},'Select a data file',defaultPath);
if fileName == 0
    fprintf(1,'Cancel button pressed. See you next time!\n')
    return
end
% Validating file selection
[~,~,fileExt] = fileparts(fileName);
switch fileExt
    case {'.smr','.smrx'}
        fprintf(1,'\nSpike2 format file\nImporting...\n')
        fStruct = dir(fullfile(filePath,strrep(fileName,fileExt,'.mat')));
        if isempty(fStruct)
            SONXimport(fopen(fullfile(filePath,fileName)))
            fprintf(1,'Done!\n')
        else
            fprintf(1,'No need to import, a mat file exists already\n')
        end
        fileName = strrep(fileName,fileExt,'.mat');
    case '.mat'
        fprintf(1,'\nMatlab file selected.\n')
    otherwise
        fprintf(1,'\nFile type not recognized. Please provide a valide file\n')
        return
end
%% Loading the variables to extract all triggers
% Loading the file into a structure
varsInFile = load(fullfile(filePath,fileName));
fNames = fieldnames(varsInFile);
chanIdx = startsWith(fNames,'chan');
headIdx = startsWith(fNames,'head');
headSub = find(headIdx);
fsArray = zeros(1,nnz(headIdx));
for chead = 1:nnz(headIdx)
    fsArray(chead) = 1e6/(varsInFile.(fNames{headSub(chead)}).sampleinterval);
    newField = erase(varsInFile.(fNames{headSub(chead)}).title,' ');
    data.(newField) = varsInFile.(strrep(fNames{headSub(chead)},'head','chan'));
end
Ns = mode(structfun(@length,data));
fs = mode(fsArray);
fprintf(1,'Sampling frequency: %.2f Hz\n',fs)

IDsignal = fieldnames(data);
% Prompting for the triggers (light/laser and piezo pulses)
idx = false(numel(fieldnames(data)),1);
[triggerIdx, iok] = listdlg(...
    'PromptString','Select all trigger signals:',...
    'ListString',IDsignal,...
    'SelectionMode','multiple',...
    'CancelString','Cancel',...
    'OKString','OK',...
    'Name','Triggers',...
    'ListSize',[160,15*numel(IDsignal)]);
if ~iok
    disp('You can always checking the file in Spike2. See you next time!')
    clearvars
    return
end
idx(triggerIdx) = true;

% Prompting for the signal(s)
[signalIdx, iok] = listdlg(...
    'PromptString','Select all signals:',...
    'ListString',IDsignal,...
    'SelectionMode','multiple',...
    'CancelString','Cancel',...
    'OKString','OK',...
    'Name','Triggers',...
    'ListSize',[160,15*numel(IDsignal(~idx))]);
if ~iok
    disp('You can always checking the file in Spike2. See you next time!')
    clearvars
    return
end
IDcontinuous = IDsignal;
IDcontinuous = IDcontinuous(signalIdx);
for ccd = 1:numel(IDcontinuous)
    ContinuousData.(IDcontinuous{ccd}) = data.(IDcontinuous{ccd});
end

% Extracting the timepoints of the step functions
cellTriggers = cell(1,numel(triggerIdx));
tSub = 1;
for cts = triggerIdx
    auxWf = StepWaveform(data.(IDsignal{cts}),fs,'Volts',IDsignal{cts});
    auxArray = auxWf.Triggers;
    if isempty(auxArray)
        % If the triggers are empty, the loop jumps to the next selected
        % trigger.
        cellTriggers(tSub) = [];
        fprintf(1,'%s skipped!\n',IDsignal{cts})
        idx(cts) = false;
        triggerIdx(tSub) = [];
        continue
    end
    Triggers.(IDsignal{cts}) = StepWaveform.SCBrownMotion(auxArray);
    cellTriggers(tSub) = {auxArray};
    tSub = tSub + 1;
end

IdxTriggers = cell2struct(cellTriggers,IDsignal(triggerIdx),2);

%% Condition search and construction
% We know that the piezo stimulation is going upwards and downwards, and
% the laser has delays with respect to the piezo onset and sometimes
% frequency stimulus.

% Piezo: Up and down and first pulse of a pulse train
% condFig = figure('Visible','off','Color',[1,1,1]);

upIdx = IdxTriggers.piezo(:,1).*data.piezo > 0;
upSub = find(upIdx);
upFst = StepWaveform.firstOfTrain(upSub/fs);
dwIdx = IdxTriggers.piezo(:,2).*data.piezo < 0;
dwSub = find(dwIdx);
dwFst = StepWaveform.firstOfTrain(dwSub/fs);

Conditions = struct('name','piezoTrainUp','Triggers',upSub(upFst));
Conditions(2).name = 'piezoTrainDw';
Conditions(2).Triggers = dwSub(dwFst);
Conditions(3).name = 'piezoTrainUpandDw';
Conditions(3).Triggers = union(upSub(upFst),dwSub(dwFst));
Conditions(4).name = 'piezoAllUp';Conditions(4).Triggers = upSub;
Conditions(5).name = 'piezoAllDw';Conditions(5).Triggers = dwSub;
Conditions(6).name = 'piezoAll';Conditions(6).Triggers = union(upSub,dwSub);

if isfield(IdxTriggers,'laser')
    disp('Development...')
    %     axes(condFig,'NextPlot','add');
    lsSub = find(IdxTriggers.laser(:,1));
    lsFst = StepWaveform.firstOfTrain(lsSub/fs, 5 - 1e-3);
    Conditions(7).name = 'laserTrain'; Conditions(7).Triggers = lsSub(lsFst);
    Conditions(8).name = 'laserAll'; Conditions(8).Triggers = lsSub;
    if isempty(Conditions(7).Triggers)
        Conditions(7) = [];
    end
    timeDelay = abs(lsSub - Conditions(3).Triggers)/fs;
    delays = uniquetol(timeDelay,0.01);
    Nc = numel(Conditions);
    for ccond = numel(Conditions)+1:numel(Conditions)+numel(delays)
        Conditions(ccond).name = sprintf('laserDelay_%d_ms',...
            floor(k*delays(ccond-...
            Nc)));
        Conditions(ccond).Triggers = lsSub(ismembertol(timeDelay,...
            delays(ccond - Nc),0.01));
    end
end

%% Spike finding
spikeFindObj = UMSDataLoader(data.Spikes,fs);
spikeFindObj.UMS2kPipeline;
spkTms = spikeFindObj.SpikeTimes;
UMSSpikeStruct = spikeFindObj.SpikeUMSStruct;
spikeFindingData = struct('thresh',[],...
    'minISI',UMSSpikeStruct.params.shadow,'spikes',spkTms,'ppms',fs/k,...
    'timestamp',[]);
if iscell(spkTms)
    spT = cell(1,numel(spkTms));
    for cn = 1:numel(spkTms)
        spT(cn) = {StepWaveform.subs2idx(round(spkTms{cn}*fs),Ns)};
    end
else
    spT = {StepWaveform.subs2idx(round(spkTms*fs), Ns)};
end
%%
% Time window surrounding the trigger [time before, time after] in seconds
timeLapse = [2.5, 5.5];

timeLapse = repmat(timeLapse,numel(Conditions),1);
timeLapse(4:6,:) = repmat([250,350]*m,3,1);
binSz = 1*m; % milliseconds or seconds
%[~,cSt] =...
%    getStacks(false(1,cols), ltOn, 'on', timeLapse * m, fs ,fs ,[],dataCell);
ERASE_kIDX = false;
IDe = [repmat('Unit ',numel(spT),1),num2str((1:numel(spT))')];
for ccon = 1:numel(Conditions)
    cn = 1;
    othNeu = setdiff(1:numel(spT),cn);
    [dSck, cSck] = getStacks(spT{cn},Conditions(ccon).Triggers,...
        'on',timeLapse(ccon,:),fs,fs, spT(othNeu),...
        struct2cell(ContinuousData));
    Na = size(dSck,3);
    kIdx = false(1,Na);
    
    [PSTH, trig, sweeps, tx_psth] =...
        getPSTH(dSck, timeLapse(ccon,:), kIdx, binSz, fs);
    [relativeSpikeTimes, tx_raster] =...
        getRasterFromStack(dSck, kIdx, false, timeLapse(ccon,:),...
        fs, ERASE_kIDX);
    expName = erase(filePath,getParentDir(filePath,1));
    plotPSTH(trig, PSTH, sweeps, binSz, timeLapse(ccon,:), expName,...
        IDe);
    plotRaster(relativeSpikeTimes, timeLapse, fs,...
        Conditions(ccon).name,erase(filePath,getParentDir(filePath,1)));
    
end
% Creation of the discrete and the continuous stacks:
%                                        First spikeing times  trigger  timeW     FS  FSt  N spiking times          continuous data.
%%
[discreteStack, contStack] = getStacks(spikeFindingData.spikes,upT,'on',timeLapse,fs,fs,spikeFindingData2.spikes,EEG.data);
%                                                               indx                   cell of indexes | INDX,  signal1, signal2,...
[Ne, Nt, Na] = size(discreteStack);
kIdx = false(1,Na);
binSz = 10; % milliseconds or seconds
ERASE_kIDX = false;
[PSTH, trig, sweeps, tx_psth] = getPSTH(discreteStack, timeLapse, kIdx, binSz, fs);
[triggeredAverageSignals, signalVariation, tx_tav] = getTriggeredAverage(contStack, kIdx, timeLapse);
[relativeSpikeTimes, tx_raster] = getRasterFromStack(discreteStack, kIdx, timeLapse, fs, ERASE_kIDX);

figure;plot(tx_psth,PSTH(1,:),tx_psth,PSTH(2,:))
plotRaster(relativeSpikeTimes, timeLapse, fs, '130626_p3');
plotRasterAsText(relativeSpikeTimes, timeLapse, fs, '130626_p3')