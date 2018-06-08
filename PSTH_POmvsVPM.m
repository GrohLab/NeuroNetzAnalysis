%% Initialization
clearvars
ToniDir = '.\Database\EphysData\AnalysisMatFiles';
expFiles = dir(fullfile(ToniDir,'*analysis.mat'));
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
%% Running through the *analysis.mat files in Toni's Directory
for cf = 1:numel(expFiles)
    [~,baseName,~]=fileparts(expFiles(cf).name);
    baseName = baseName(1:end-8);
    if mkdir(ToniDir,baseName)
        if ~exist(fullfile(ToniDir,baseName,expFiles(cf).name),'file')
            disp(['Trying to move file ',expFiles(cf).name])
            if ~movefile(fullfile(ToniDir,expFiles(cf).name),...
                    fullfile(ToniDir,basename))
                writeLogFile(logFile,[expFiles(cf).name, ' couldn''t be found'])
                continue;
            end
        end
        fileName = fullfile(ToniDir,basename,expFiles(cf).name);
        try
            load(fileName,'filteredResponse','Triggers','EEG')
            if ~isempty(filteredResponse.data)
                spkWvf = SpikeWaveform(filteredResponse.data,fs);
                spkWvf.getSpikes_UMS;
                spT = spkWvf.Triggers;
                UMSSpikeStruct = spkWvf.UMSdata.SpikeUMSStruct;
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
        try
            if UMS
                save(fileName,'UMSSpikeStruct','a')
                writeLogFile(logFile,[fileName,...
                ' UMS2k structure appended'])
            end
        catch
            writeLogFile(logFile,[fileName, ' unwritable!'])
            disp([fileName,' unwritable!'])
        end
        cellType = RecDB(baseName,'PhysioNucleus');
        cellType = string(cellType.PhysioNucleus);
        timeLapse = [500,1e3]*m;
        % Whisker onset
        whiskWf = StepWaveform(Triggers.whisker,fs);
        RaF = whiskWf.Triggers;
        whiskerMovement = Triggers.whisking;
        fieldNames = fieldnames(Triggers);
        allEvents = struct2cell(Triggers);
        consEvents = allEvents(~ismember(fieldNames,{'whisking','whisker'}));
        [PSTHstackW,LFPstackW,WstackW] = getStack(spT,RaF,'on',...
            timeLapse,fs,EEG.data,Triggers.whisking,consEvents);
        [~,FigID] = plotTriggeredEvents(PSTHstackW,LFPstackW,WstackW,timeLapse,cellType,0.01,fs);
        
        % Whisker offset
        [PSTHstackNW,LFPstackNW,WstackNW] = getStack(spT,RaF,'off',...
            timeLapse,fs,EEG.data,Triggers.whisking,consEvents);
        

    else
        writeLogFile(logFile,[fullfile(ToniDir,baseName), ' couldn''t be created'])
        disp([fullfile(ToniDir,baseName), ' couldn''t be created'])
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