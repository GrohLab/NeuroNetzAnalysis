%% Initialization
clearvars
ToniDir = '.\Database\EphysData\AnalysisMatFiles';
expFiles = dir(fullfile(ToniDir,'*_00*'));
UMS = false;
logFile = 'TONI-LOG.txt';
m = 1e-3;
fs = 2e4;

EphysPath = '.\Database\EphysData\';
DBPath = getParentDir(EphysPath,1);
try
    load([DBPath, 'ephys_database.mat'])
catch ME
    disp(ME.getReport)
    disp('The Database couldn''t be loaded!')
end
RecDB.Properties.RowNames = strcat(RecDB.AnimalName,{'_'},...
    datestr(RecDB.StartTime, 'HH_MM_SS'));
try
    load([EphysPath,'DiscreteData.mat'], 'DiscreteData');
catch ME
    disp(ME.getReport)
    disp('The discrete data couldn''t be loaded')
end
%% For Loop body
for cf = 1:numel(expFiles)
    %% File managing section
    % Running through the folders or files in the given path and 'organize'
    % the files in homonime folders. Then the algorithm saves the results
    % together with the data file.
    if expFiles(cf).isdir
        analysisFile = dir(fullfile(ToniDir,expFiles(cf).name,'*analysis.mat'));
        analysisFile = fullfile(analysisFile.folder,analysisFile.name);
        [~,baseName,~] = fileparts(analysisFile);
    else
        [~,baseName,~]=fileparts(expFiles(cf).name);
        analysisFile = fullfile(ToniDir,expFiles(cf).name);
    end
    baseName = baseName(1:end-8);
    if mkdir(ToniDir,baseName)
        expDir = fullfile(ToniDir,baseName);
        if ~exist(analysisFile,'file')
            disp(['Trying to move file ',expFiles(cf).name])
            if ~movefile(fullfile(ToniDir,expFiles(cf).name),...
                    fullfile(ToniDir,baseName))
                writeLogFile(logFile,[expFiles(cf).name, ' couldn''t be found'])
                continue;
            end
        end
        fileName = analysisFile;
        try
            % Load the spike times from the UMS structure and use it
            % directly to compute the stack.
            if exist('UMSSpikeStruct','var')
                clearvars UMSSpikeStruct
            end
            load(fileName,'UMSSpikeStruct') % Removed a 't' to create an error
            % spT = UMSSpikeStruct.spiketimes; % Indexes
            UMSObject = UMSDataLoader();
            UMSObject.changeUMSStructure(UMSSpikeStruct)
            UMSObject.getSpikeTimes
            spT = UMSObject.SpikeTimes; % Indexes
            load(fileName,'Triggers','EEG')
            UMS = false;
        catch
            try
                load(fileName,'filteredResponse','Triggers','EEG')
                if ~isempty(filteredResponse.data)
                    spkWvf = SpikeWaveform(filteredResponse.data,fs);
                    spkWvf.getSpikes_UMS;
                    spT = spkWvf.Triggers;
                    UMSSpikeStruct = spkWvf.UMSdata.SpikeUMSStruct;
                    h = spkWvf.plot;
                    try
                        while strcmp(h.BeingDeleted,'off')
                            waitforbuttonpress
                        end
                    catch
                    end
                    UMS = true;
                elseif ~isempty(UMS)
                    writeLogFile(logFile,[baseName,...
                        ' no filtered response found.\t Using Toni''s extracted time stamps'])
                    expIdx = find(strcmp(RecDB.Properties.RowNames, baseName));
                    spT = DiscreteData(expIdx).Spikes;
                    UMS = false;
                end
            catch
                writeLogFile(logFile,[expFiles(cf).name,...
                    ' unreadable!'])
                continue
            end
        end
        try
            if UMS
                save(fileName,'UMSSpikeStruct','-append')
                writeLogFile(logFile,[fileName,...
                ' UMS2k structure appended'])
            end
        catch
            writeLogFile(logFile,[fileName, ' unwritable!'])
            disp([fileName,' unwritable!'])
        end
        %% Stack creation for whisker on and off-set
        % Initialization for the stack building
        cellType = RecDB(baseName,'PhysioNucleus');
        cellType = string(cellType.PhysioNucleus);
        timeLapse = [100,300]*m;
        whiskWf = StepWaveform(Triggers.whisker,fs);
        RaF = whiskWf.Triggers;
        whiskerMovement = Triggers.whisking;
        fieldNames = fieldnames(Triggers);
        allEvents = struct2cell(Triggers);
        consEvents = allEvents(~ismember(fieldNames,{'whisking','whisker'}));
        condNames = fieldNames(~ismember(fieldNames,{'whisking','whisker'}));
        ConditionsTest = cell(1,numel(consEvents));
        for ccon = 1:numel(consEvents)
            auxStpWvf = StepWaveform(consEvents{ccon},fs,'mV',condNames{ccon});
            condTriggers = auxStpWvf.Triggers;
            auxStruct = struct('name',condNames{ccon},...
                'Triggers',condTriggers,...
                'delay',0);
            ConditionsTest(ccon) = {auxStruct};
        end
        binningTime = 10*m;
        fsLFP = 1e3;
        % Overwritting the conditions test variable with each run. This is
        % time consuming, so take that into consideration.
        save(fileName,'ConditionsTest','-append')
        fprintf('File %s saved!\n',fileName)
        % Whisker onset
        [PSTHstackW,contStackW] = getStacks(spT, RaF, 'on',...
            timeLapse, fs,fsLFP,consEvents, EEG.data, Triggers.whisking);
%         [~,FigIDW,kO] =...
%             plotTriggeredEvents(PSTHstackW, squeeze(contStackW(1,:,:)),...
%             squeeze(contStackW(2,:,:)),...
%             timeLapse, cellType, binningTime, fs,fsLFP);
        [~,FigIDW,kO] =...
            plotTriggeredEvents(PSTHstackW, (contStackW(1,:,:)),...
            (contStackW(2,:,:)),...
            timeLapse, cellType, binningTime, fs,fsLFP);
        savePDF_FIG(FigIDW,'_W',expDir,baseName)
%         if ~exist(fullfile(expDir,[baseName,'_W.pdf']),'file')
%             print(FigIDW,fullfile(expDir,[baseName,'_W.pdf']),...
%                 '-dpdf','-fillpage')
%             savefig(FigIDW,fullfile(expDir,[baseName,'_W.fig']))
%         end
        close(FigIDW)
        % Whisker offset
        [PSTHstackNW,contStackNW] = getStacks(spT, RaF, 'off',...
            timeLapse, fs,fsLFP,consEvents, EEG.data, Triggers.whisking);
        % [~,FigIDNW] =...
        % plotTriggeredEvents(PSTHstackNW, squeeze(contStackNW(1,:,:)),...
        % squeeze(contStackNW(2,:,:)),...
        % timeLapse, cellType, binningTime, fs,fsLFP);
        [~,FigIDNW] =...
            plotTriggeredEvents(PSTHstackNW, (contStackNW(1,:,:)),...
            (contStackNW(2,:,:)),...
            timeLapse, cellType, binningTime, fs,fsLFP);
        savePDF_FIG(FigIDNW,'_NW',expDir,baseName)
%         if ~exist(fullfile(expDir,[baseName,'_NW.pdf']),'file')
%             print(FigIDNW,fullfile(expDir,[baseName,'_NW.pdf']),...
%                 '-dpdf','-fillpage')
%             savefig(FigIDNW,fullfile(expDir,[baseName,'_NW.fig']))
%         end
        close(FigIDNW)
        % Sanity plot
        figIDS = figure('Name',[char(cellType), ' considered events'],...
            'Color',[1,1,1]);
        try
            subplot(2,1,1);mesh(squeeze(sum(PSTHstackW(3:end,:,:),2)));ylabel('Event')
            box off;xlabel('Alignment point');zlabel('Count');
            title('Segmented events in stack')
            subplot(2,1,2);mesh(squeeze(sum(PSTHstackW(3:end,:,~kO),2)))
            title('Events to include in the PSTH');ylabel('Event');box off
            xlabel('Alignment point');zlabel('Count');
        catch
            disp('Maybe there is no ''clean'' data')
            writeLogFile(logFile,[fullfile(ToniDir,baseName), ' sanity plot not created'])
        end
        savePDF_FIG(figIDS,'_consideredEvents',expDir,baseName)
%         if exist(fullfile(expDir,[baseName,'consideredEvents.pdf']),'file')
%             print(figIDS,fullfile(expDir,[baseName,'consideredEvents.pdf']),...
%                 '-dpdf','-fillpage')
%             savefig(figIDS,fullfile(expDir,[baseName,'consideredEvents.fig']))
%         end
        close(figIDS)
    else
        writeLogFile(logFile,[fullfile(ToniDir,baseName), ' couldn''t be created'])
        disp([fullfile(ToniDir,baseName), ' couldn''t be created'])
    end
end

function savePDF_FIG(h,postfix,path,bName)
if ~exist(fullfile(path,[bName,postfix,'.pdf']),'file') ||...
        ~exist(fullfile(path,[bName,postfix,'.fig']),'file')
    print(h,fullfile(path,[bName,postfix,'.pdf']),...
        '-dpdf','-fillpage')
    savefig(h,fullfile(path,[bName,postfix,'.fig']))
end
end

function myStatus = writeLogFile(logFileName,myMessage)
fid = fopen(logFileName,'a');
if fid == -1
    disp('----------Log file unavailable/corrupt!------------')
    myStatus = false;
else
    fprintf(fid,'%s: %s\n', datestr(now, 0), myMessage);
    myStatus = true;
end
end