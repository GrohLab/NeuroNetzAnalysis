%% Constant values
% Prefixes
anaFlag = false;
k = 1e3; % Kilo
m = 1e-3; % milli
zs = @(x)mean(x)/std(x);
addFst = @(x,y)cat(find(size(x)~=1),y,x);
addLst = @(x,y)cat(find(size(x)~=1),x,y);
fstCell = @(x) x{1};
invalidCharacters = {'+','-'};
uninterestingVals = @(x)~(isempty(x) || isnan(x) || isinf(x) || x == 0);
findUnpairedPulse = @(x) cat(1,false,reshape(diff(sort(x)) == 0,numel(x)-1,1));
printFig = @(x,fname) set(x,'RendererMode','manual','Renderer','painters',...
    'PaperOrientation','landscape','Name',fname);
%% Select data file
defaultPath = 'E:\Data\Jittering'; % Path where you would like to open a windows explorer window to select your data file
% defaultPath = 'C:\Users\jefe_\Documents\SampleData';
[fileName, filePath] =...
    uigetfile({'*analysis.mat';'*.mat';'*.smr';'*.smrx'},...
    'Select a data file',defaultPath);
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
    case {'.mat','analysis.mat'}
        fprintf(1,'\nMatlab file selected.\n')
        anaFlag = any(strfind(fileName,'analysis.mat'));
        if anaFlag
            fprintf(1,'Skipping to the plotting section... \n')
        end
    otherwise
        fprintf(1,'\nFile type not recognized. Please provide a valide file\n')
        return
end

if ~anaFlag
    [~,~,fExt] = fileparts(fileName);
    afn = fullfile(filePath,[erase(fileName,fExt),'_analysis.mat']);
end
%% Extracting the information from the selected file
if ~anaFlag
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
        newField = erase(varsInFile.(fNames{headSub(chead)}).title,{' '});
        if contains(newField,invalidCharacters)
            newFieldProposal = strrep(newField,invalidCharacters,{'NVS'});
            wrongName = strfind(newFieldProposal,'NVS');
            newField = newFieldProposal{~cellfun(@isempty,wrongName)};
        end
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
        'ListSize',[160,15*numel(IDsignal)]);
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
        Triggers.(IDsignal{cts}) =...
            StepWaveform.SCBrownMotion(auxArray) ~= 0;
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
    
    upIdx = IdxTriggers.piezo .* repmat(data.piezo,1,2) > 0;
    dwIdx = IdxTriggers.piezo .* repmat(data.piezo,1,2) < 0;
    pzIdx = upIdx | flip(dwIdx,2);
    pzSub = [find(pzIdx(:,1)),find(pzIdx(:,2))];
    pzFst = StepWaveform.firstOfTrain(pzSub(:,1)/fs);

    upSub = [find(upIdx(:,1)),find(upIdx(:,2))];
    upFst = StepWaveform.firstOfTrain(upSub(:,1)/fs);
    dwSub = [find(dwIdx(:,2)), find(dwIdx(:,1))];
    dwFst = StepWaveform.firstOfTrain(dwSub(:,2)/fs);
   
    Conditions = struct('name','Piezo','Triggers',union(upSub,dwSub,'rows'));
    
    if isfield(IdxTriggers,'laser')
        % Searching for frequency stimulation
        lsSub = find(IdxTriggers.laser(:,1));
        lsIpi = diff(lsSub/fs);
        lsFst = StepWaveform.firstOfTrain(lsSub/fs, 5 - 1e-3);
        
        lsIpi = lsIpi(~lsFst(1:end-1));
        pulsFreq = 1./diff(lsSub./fs);
        freqCond = round(uniquetol(pulsFreq,0.1/max(pulsFreq)));
        freqCond = unique(freqCond(freqCond > 0));
        fprintf(1,'Frequency stimulation:')
        
        if isempty(freqCond) || numel(freqCond) == 1
            freqCond = 0;
            fprintf(' None')
        elseif numel(freqCond) > 1
            Nfre = numel(freqCond);
            pulsFreq = 1./diff(lsSub./fs);
            lsuSub = lsSub(lsFst);
            lsdSub = find(IdxTriggers.laser(:,2));
            lsLst = flip(StepWaveform.firstOfTrain(flip(lsdSub)/fs, 5-1e-3));
            lsdSub = lsdSub(lsLst);
            lsCon = [lsuSub,lsdSub];
            lsIdx = false(size(lsCon,1),Nfre);
            for cdl = 1:Nfre
                fprintf(1,' %.2f',freqCond(cdl))
            end
        end
        fprintf(1,' Hz\n')
        
        % Searching for delays in the data with respect to the piezo
        
        maxPulses = min(numel(lsSub),size(Conditions.Triggers,1));
        dm = distmatrix(lsSub/fs,Conditions.Triggers(:,1)/fs);
        [srtDelay, whr] = sort(dm(:),'ascend');
        [lghtSub, piezSub] = ind2sub(size(dm),whr(1:maxPulses));
        timeDelay = srtDelay(1:maxPulses);
        delays = 10.^uniquetol(log10(timeDelay),0.01/log10(max(abs(timeDelay))));
        delays(delays > 1) = [];
        if std(delays.*1e3) < 1
            delays = mean(delays);
        end
        Ndel = numel(delays);
        fprintf(1,'Delays found:')
        lsDel = false(length(lsSub),Ndel);
        Ncond = numel(Conditions);
        for cdl = 1:Ndel
            fprintf(1,' %.1f',delays(cdl)*1e3)
            lsDel(:,cdl) = ismembertol(log10(timeDelay),log10(delays(cdl)),...
                abs(0.01/log10(max(delays))));
            Conditions(Ncond + cdl).name = sprintf('Delay %0.3f s',...
                delays(cdl));
            Conditions(Ncond + cdl).Triggers =...
                Conditions(1).Triggers(sort(piezSub(lsDel(:,cdl))),1); 
        end
        fprintf(1,' ms\n')
        [~,lghtSub] = min(dm,[],2,'omitnan');
        [~,piezSub] = min(dm,[],1,'omitnan');
        piezSub = piezSub';
        loneLaser = findUnpairedPulse(lghtSub);
        lonePiezo = reshape(diff(sort(piezSub)) == 0,numel(piezSub)-1,1);
        Conditions(Ncond + cdl + 1).name = 'Laser Control';
        Conditions(Ncond + cdl + 1).Triggers = lsOn(find(loneLaser));
        Conditions(Ncond + cdl + 2).name = 'Puff Control';
        Conditions(Ncond + cdl + 2).Triggers = pzUp(find(lonePiezo));
    end
    
    %% Spike finding
    spikeFindObj = UMSDataLoader(data.Spikes,fs);
    spikeFindObj.UMS2kPipeline;
    spkTms = spikeFindObj.SpikeTimes;
    UMSSpikeStruct = spikeFindObj.SpikeUMSStruct;
    spikeFindingData = struct('thresh',[],...
        'minISI',UMSSpikeStruct.params.shadow,'spikes',spkTms,'ppms',fs/k,...
        'timestamp',[]);
    %% Save *analysis.mat file
    
    Head = struct('FileInfo',fullfile(filePath,fileName),...
        'IDs',{IDsignal},'tSub',triggerIdx,'sSub',signalIdx,'Samples',Ns,...
        'SamplingFrequency',fs);
    save(afn,'Triggers','ContinuousData','Conditions','UMSSpikeStruct','Head')
    
else
    varsInFile = load(fullfile(filePath,fileName));
    
    Ns = varsInFile.Head.Samples;
    fs = varsInFile.Head.SamplingFrequency;
    IDsignal = varsInFile.Head.IDs;
    triggerIdx = varsInFile.Head.tSub;
    signalIdx = varsInFile.Head.sSub;
    
    Triggers = varsInFile.Triggers;
    Conditions = varsInFile.Conditions;
    ContinuousData = varsInFile.ContinuousData;
    UMSSpikeStruct = varsInFile.UMSSpikeStruct;
    
    clearvars varsInFile
    
    UMSobj = UMSDataLoader();
    UMSobj.changeUMSStructure(UMSSpikeStruct);
    spkTms = UMSobj.SpikeTimes;
end
%%
if iscell(spkTms)
    spT = cell(numel(spkTms),1);
    for cn = 1:numel(spkTms)
        spT(cn) = {StepWaveform.subs2idx(round(spkTms{cn}*fs)',Ns)};
    end
else
    spT = {StepWaveform.subs2idx(round(spkTms*fs)', Ns)};
end

%%
% Time window surrounding the trigger [time before, time after] and bin
% size for the PSTH. Both quantities in seconds.
timeLapse = [350, 350]*m;
binSz = 10*m;

ERASE_kIDX = false;
clID = UMSDataLoader.getClustersID(UMSSpikeStruct,'good');
subOffst = numel(clID) - 1;
IDe = IDsignal(...
    ~cellfun(@strcmpi,IDsignal,repmat({'Spikes'},numel(IDsignal),1)));
rID = num2cell([repmat('Cluster ',numel(clID),1),num2str(clID)],2);


ntSub = triggerIdx + subOffst;
nsSub = signalIdx + subOffst;
consideredSignalsIdx = false(size(IDe));

consEvnts = struct2cell(Triggers);
consEvnts(tIdx) = [];
eID = fieldnames(Triggers);
tIdx = strcmp(eID,{'piezo'});
othNeu = setdiff(1:numel(spT),1);

if ~(isempty(spT(othNeu)) || isempty(consEvnts))
    consEvnts = cat(1,spT(othNeu),consEvnts);
end
OVW = false; % Overwrite figure flag

expName = [erase(erase(filePath,getParentDir(filePath,1)),filesep),'->',...
        erase(fileName,'.mat')];

trigSubs = pzSub;
if contains(fileName,'long')
    trigSubs = pzSub(pzFst,:);
else
    disp('Not long experiment...')
end
[dStack, cStack] = getStacks(spT{1},trigSubs,'on',timeLapse,fs,fs,...
    consEvnts,struct2cell(ContinuousData));

[Ne,~,Na] = size(dStack);

kIdx = false(1, Na);
interestingEvents = true(1,Ne-2);
%%
eIDx = purgeTrials(dStack,timeLapse,tIdx,~tIdx,...
    false(size(tIdx)),eID,fs);
if ~sum(eIDx)
    eIDx = true(Na,1);
end
[PSTH, trig, sweeps] =...
    getPSTH(dStack, timeLapse, ~eIDx, binSz, fs);
koIdx = true(size(PSTH,1),1);
koIdx(tIdx) = false;
prID = rID;
prID(1) = [];
plotPSTH(trig, PSTH, sweeps, binSz, timeLapse, expName, [prID;eID(~tIdx)], koIdx,...
    padarray(tIdx,size(koIdx,1)-size(tIdx,1),'post'), fs);

[relativeSpikeTimes] =...
    getRasterFromStack(dStack, ~eIDx, interestingEvents, timeLapse,...
    fs, ERASE_kIDX);
plotRaster(relativeSpikeTimes,timeLapse,fs,['Raster for ',expName],[rID;eID(~tIdx)]);
    
%%

if anaFlag
    endName = erase(fileName,'analysis.mat');
else
    endName = erase(fileName,'.mat');
end
for ccon = 1:numel(Conditions)
    [dStack, cSck] = getStacks(...
        spT{1},... Main neuron
        Conditions(ccon).Triggers, 'on',... Triggers, onset
        timeLapse(ccon,:), fs, fs,... Time lapse, and sampling frequencies
        consEvnts,... Other discrete events in the experiment.
        struct2cell(ContinuousData)); % Continuous data
    
    
    
    
    
    [PSTH, trig, sweeps] =...
        getPSTH(dStack, timeLapse(ccon,:), kIdx, binSz(ccon), fs);
    [relativeSpikeTimes] =...
        getRasterFromStack(dStack, kIdx, interestingEvents, timeLapse(ccon,:),...
        fs, ERASE_kIDX);
    expName = [erase(erase(filePath,getParentDir(filePath,1)),filesep),'->',...
        endName];

    psthFig = plotEasyPSTH(trig, PSTH, sweeps, binSz(ccon),...
        timeLapse(ccon,:), fs);
    set(psthFig.Children.Title,'String',[expName,' ',Conditions(ccon).name])
    printFig(psthFig,sprintf('PSTH for %s',expName))
    psthFigName = fullfile(filePath,[expName,'_PSTH_',Conditions(ccon).name,'.pdf']);
    if ~exist(psthFigName,'file') || OVW
        print(psthFig,psthFigName,'-dpdf','-bestfit')
    end
    rastFig = plotRaster(relativeSpikeTimes, timeLapse(ccon,:), fs,...
        [expName,' ',Conditions(ccon).name],rID);
    printFig(rastFig,sprintf('Raster for %s',expName))
    rastFigName = fullfile(filePath,[expName,'_Raster_',Conditions(ccon).name,'.pdf']);
    if ~exist(rastFigName,'file') || OVW
        print(rastFig, rastFigName,'-dpdf','-bestfit')
    end
end