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
        lsSub = find(IdxTriggers.laser(:,1));
        lsIpi = diff(lsSub/fs);
        lsFst = StepWaveform.firstOfTrain(lsSub/fs, 5 - 1e-3);
        
        lsIpi = lsIpi(~lsFst(1:end-1));
        freqCond = round(uniquetol(1./lsIpi,0.01));
        freqCond = freqCond(freqCond > 0);
        fprintf(1,'Frequency stimulation:')
        if isempty(freqCond)
            freqCond = 0;
        else
            for cdl = 1:numel(freqCond)
                fprintf(1,' %.2f',freqCond(cdl))
            end
            fprintf(1,'\n')
        end
        try
            if isempty(lsFst) || ~sum(lsFst)
                timeDelay = abs(lsSub - Conditions(3).Triggers)/fs;
            else
                timeDelay = abs(lsSub(lsFst) - Conditions(3).Triggers)/fs;
            end
        catch TDE
            fprintf(1,'Seems that a pulse was truncated...\n')
            fprintf(1,'----Please pay attention when you finish recording!\n')
            if isempty(lsFst) || ~sum(lsFst)
                maxPulses = min(numel(lsSub),numel(Conditions(3).Triggers));
                timeDelay = abs(lsSub(1:maxPulses) -...
                    Conditions(3).Triggers(1:maxPulses))/fs;
            else
                maxPulses = min(sum(lsFst),numel(Conditions(3).Triggers));
                dm = distmatrix(lsSub(lsFst)/fs,Conditions(3).Triggers/fs);
                srtDelay = sort(dm(:),'ascend');
                timeDelay = srtDelay(1:maxPulses);
            end
        end
        delays = uniquetol(timeDelay,0.01);
        if std(delays.*1e3) < 1
            delays = mean(delays);
        end
        fprintf(1,'Delays found:')
        for cdl = 1:numel(delays)
            fprintf(1,' %.2f',delays(cdl))
        end
        fprintf(1,'\n')
        Nc = numel(Conditions);
        Ndel = numel(delays);
        Nfre = numel(freqCond);
        for ccond = Nc+1:Nc*(Ndel*Nfre+1)
            subFreq = ceil(((ccond/Nc)-1)*(1/Ndel));
            subDel = ceil(((ccond/Nc)-1)*(1/Nfre));
            subCond = mod(ccond-1,Nc)+1;
            Conditions(ccond).name = cat(2,Conditions(subCond).name,...
                sprintf('+laser(%d ms & %d Hz)',round(k*delays(subDel)),...
                freqCond(subFreq)));
            % TODO: Find the relevant laser + piezo combination!
            dm = log(distmatrix(lsSub/fs,Conditions(subCond).Triggers/fs));
            [y,x] = find(ismembertol(dm,min(dm(:)),0.001));
            Conditions(ccond).Triggers = Conditions(subCond).Triggers(x);
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
timeLapse = [0, 50*m];
binSz = 50*m; 

%[~,cSt] =...
%    getStacks(false(1,cols), ltOn, 'on', timeLapse * m, fs ,fs ,[],dataCell);
ERASE_kIDX = false;
clID = UMSDataLoader.getClustersID(UMSSpikeStruct,'good');
subOffst = numel(clID) - 1;
IDe = IDsignal(...
    ~cellfun(@strcmpi,IDsignal,repmat({'Spikes'},numel(IDsignal),1)));
IDe = [num2cell([repmat('Cluster ',numel(clID),1),num2str(clID)],2);IDe];
rID = [num2cell([repmat('Cluster ',numel(clID),1),num2str(clID)],2);...
    fieldnames(Triggers)];
ntSub = triggerIdx + subOffst;
nsSub = signalIdx + subOffst;
consideredSignalsIdx = false(size(IDe));

othNeu = setdiff(1:numel(spT),1);
if isempty(spT(othNeu))
    consEvnts = struct2cell(Triggers);
else
    consEvnts = cat(1,spT(othNeu),struct2cell(Triggers));
end
OVW = false; % Overwrite figure flag
%%
for ccon = 1:numel(Conditions)
    [dSck, cSck] = getStacks(...
        spT{1},... Main neuron
        Conditions(ccon).Triggers, 'on',... Triggers, onset
        timeLapse(ccon,:), fs, fs,... Time lapse, and sampling frequencies
        consEvnts,... Other discrete events in the experiment.
        struct2cell(ContinuousData)); % Continuous data
    [Ne,~,Na] = size(dSck);
    
    kIdx = false(1,Na);
    interestingEvents = true(1,Ne-2);
    
    [PSTH, trig, sweeps] =...
        getPSTH(dSck, timeLapse(ccon,:), kIdx, binSz(ccon), fs);
    [relativeSpikeTimes] =...
        getRasterFromStack(dSck, kIdx, interestingEvents, timeLapse(ccon,:),...
        fs, ERASE_kIDX);
    expName = [erase(erase(filePath,getParentDir(filePath,1)),filesep),...
        erase(fileName,'.mat')];
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