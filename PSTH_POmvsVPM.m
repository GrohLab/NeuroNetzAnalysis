%% Initialization
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
                    ' no filtered response found.'])
                spT = find(strcmp(RecDB.Properties.RowNames, baseName));
                
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
            end
        catch
            writeLogFile(logFile,[fileName,...
                ' unwritable!'])
        end
        cellType = RecDB(baseName,'PhysioNucleus');
        % Whisker onset
        RaF = whiskWf.Triggers;
        whiskerMovement = Triggered
        [PSTHstack1,LFPstack1,~] = getPSTH(spT,Raf(:,1),...
            [100,200]*m,fs,EEG.data,{});
        % Whisker offset
        whiskWf = StepWaveform(Triggers.whisker,fs);
        [PSTHstack2,LFPstack2,~] = getPSTH(spT,RaF(:,2),[100,200]*m,...
            fs,EEG.data,{});
        
        
    else
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