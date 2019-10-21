function iOk = loadTriggerData(dataDir)
iOk = false;
binFiles = dir(fullfile(dataDir,'*.bin'));
[~,expName,~] = fileparts( binFiles.name);
% Loading the sampling frequency, the sorted clusters, and the conditions
% and triggers.
expSubfix = fullfile(dataDir,expName);
assignin('base','expSubfix',expSubfix)
try
    load([expSubfix,'_sampling_frequency.mat'],'fs')
catch
    fprintf(1,'Seems like the kilosort-phy pipeline hasn''t been touched!\n')
    return
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
iOk = true;
end