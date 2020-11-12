function iOk = loadTriggerData(dataDir)
iOk = false;
binFiles = dir(fullfile(dataDir,'*.bin'));
smrxFiles = dir(fullfile(dataDir,'*.smrx'));
try
[~,expName,~] = fileparts( binFiles(1).name);
catch
    fprintf(1,'No binary file in the folder!\n')
    try
       [~,expName,~] = fileparts( smrxFiles(1).name); 
    catch
        fprintf(1,'No smrx file in the folder either!\n')
        forthAns =...
            questdlg('No binary nor smrx files in this folder! Continue?',...
            'Continue?','Yes','No','No');
        if strcmpi(forthAns,'No')
            return
        end
    end
end
% Loading the sampling frequency, the sorted clusters, and the conditions
% and triggers.
expSubfix = fullfile(dataDir,expName);

% chanMap = readNPY(fullfile(dataDir,'channel_map.npy'));
% chanPos = readNPY(fullfile(dataDir,'channel_positions.npy'));

% assignin('base','chanMap',chanMap)
% assignin('base','chanPos',chanPos)

assignin('base','expSubfix',expSubfix)
assignin('base','expName',expName)
try
    load([expSubfix,'_sampling_frequency.mat'],'fs')
catch
    try
        [~,smrxBaseName] = fileparts(smrxFiles(1).name);
        load([fullfile(dataDir,smrxBaseName),'_sampling_frequency.mat'],'fs')
    catch
        fprintf(1,'Did not find the sampling frequency. Cannot continue\n')
        return
    end
end
assignin('base','fs',fs)
try
    load([expSubfix,'analysis.mat'],'Conditions','Triggers')
catch
    if exist([expSubfix, '_CondSig.mat'],'file')
        getDelayProtocol(dataDir);
    else
        try
            getConditionSignalsBF(fopen([expSubfix,'.smrx']))
            getDelayProtocol(dataDir);
        catch
            fprintf(1,'Confusing naming. Cannot continue\n')
            return
        end
    end
    load([expSubfix,'analysis.mat'],'Conditions','Triggers')
end
assignin('base','Conditions',Conditions)
assignin('base','Triggers',Triggers)
try
    load([expSubfix,'_all_channels.mat'],'sortedData')
catch
    try
        importPhyFiles(dataDir);
    catch
        fprintf(1,'Error importing the phy files into Matlab format\n')
        return
    end
    load([expSubfix,'_all_channels.mat'],'sortedData')
end
assignin('base','sortedData',sortedData)
try
    clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
catch
    fprintf(1, 'No ''cluster_info.tsv'' file found!\n')
    fprintf(1, 'Be sure to save your sorting with phy\n')
    clInfo = table();
end
assignin('base','clInfo',clInfo)
iOk = true;
end