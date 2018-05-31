ToniDir = '.\Database\EphysData\AnalysisMatFiles';
fullfile(ToniDir,'*analysis.mat')
UMS = false;
logFile = 'TONI-LOG.txt';
for cf = 1:numel(expFiles)
    [~,baseName,~]=fileparts(expFiles(cf).name);
    baseName = baseName(1:end-8);
    if mkdir(ToniDir,baseName)
        if movefile(fullfile(ToniDir,expFiles(cf).name),...
                fullfile(ToniDir,basename))
            fileName = fullfile(ToniDir,basename,expFiles(cf).name);
            try
                load(fileName,'filteredResponse','spikeFindingData',...
                    'Triggers','EEG')
                if ~isempty(filteredResponse.data)
                    spkWvf = SpikeWaveform(filteredResponse.data,2e4);
                    spkWvf.getSpikes_UMS;
                    spT = spkWvf.Triggers;
                    UMSSpikeStruct = spkWvf.UMSdata.SpikeUMSStruct;
                    UMS = true;
                elseif ~isempty(spikeFindingData.spikes)
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
            [PSTHstack,LFPstack,~] = getPSTH(spT,Tirggers.whisker,...
                [100,200]*1e-3,2e4,EEG.data,{});
        else
        end
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