CFA_DIR = 'E:\LFP_Frequency_Analysis\Juxta\CFA\';
SAL_DIR = 'E:\LFP_Frequency_Analysis\Juxta\Saline\';
%%
dirFiles = dir([CFA_DIR,'*.mat']);
CFA_LFP = cell(length(dirFiles),4);
CFA_SPCT = cell(length(dirFiles),2);
CFA_MECH = cell(length(dirFiles),1);
stmSub = [];
mechSub=[];
figure('Name','CFA power spectrum','color',[1,1,1]);
% function outputArg1 = processPainLFP(directory)
%% Folder processing
inDirs = dir(directory);
ccon = 1;
while ccon <= numel(inDirs)
    if inDirs(ccon).name == "." || inDirs(ccon).name == ".." || ...
            ~inDirs(ccon).isdir
        inDirs(ccon) = [];
        disp('File or dot(s) found')
        continue;
    end
    fprintf('Processing %s folder...\n',inDirs(ccon).name)
    lfpStack = [];
    dirFiles = dir(fullfile(inDirs(ccon).folder,inDirs(ccon).name,'*analysis.mat'));
    for cf = 1:length(dirFiles)
        %% Loading and a small preprocessing
        load(fullfile(dirFiles(cf).folder,dirFiles(cf).name),...
            'EEG','Conditions')
        fs = EEG.header.SamplingFrequency;
        desFs = 1e3;
        lfpAux = EEG.data;
        if ~isa(lfpAux,'double')
            lfpAux = double(lfpAux);
        end
        if iscolumn(lfpAux)
            lfpAux = lfpAux';
        end
        %% Resampling the original signal
        N = length(lfpAux);
        tx = (0:N-1)*(1/fs);
        desTx = 0:1/desFs:tx(end);
        lfpTS = timeseries(lfpAux,tx);
        lfpTS = lfpTS.resample(desTx);
        lfpAux = squeeze(lfpTS.Data);
        lfpAux = brainwaves(lfpAux,desFs,{'alpha',1e-2,200});
        lfpAux = iir50NotchFilter(lfpAux,desFs);
        N = length(lfpAux);
        
        %% Transforming the triggers' time points
        stmSub = [];
        for cc = 1:numel(Conditions)
            trgSub = Conditions{cc}.Triggers;
            trgSub = round(trgSub .* (fs/desFs));
            
            disp(Conditions{cc}.name)
            if isrow(trgSub)
                trgSub = trgSub';
            end
            stmSub =[stmSub; trgSub, trgSub + round(5.35*desFs)];
            if strcmp(Conditions{cc}.name,'MechControl')
                mecSub = [trgSub, trgSub+ round(5.35*desFs)];
            end
        end
        stmSub = round(stmSub);
        mecSub = round(mecSub);
        %% Get Stacks
        [~,stackAux] =...
            getStacks(false(1,length(lfpAux)),...
            mecSub,'on',[1,6],desFs,desFs,[],lfpAux);
        lfpStack = cat(3,lfpStack,zeros(1,size(stackAux,2),1),stackAux);
    end
    lfpConditionsStack(ccon) = {lfpStack};
    ccon = ccon + 1;
end

%%
dirFiles = dir([SAL_DIR,'*.mat']);
SAL_LFP = cell(length(dirFiles),4);
SAL_SPCT = cell(length(dirFiles),2);
SAL_MECH = cell(length(dirFiles),1);
figure('Name','Saline power spectrum','Color',[1,1,1]);
for cf = 1:length(dirFiles)
    load(fullfile(SAL_DIR,dirFiles(cf).name),'EEG','Conditions')
    fs = EEG.header.SamplingFrequency;
    lfpAux = brainwaves(double(EEG.data'),fs,{'alpha',1e-2,200});
    lfpAux = iir50NotchFilter(lfpAux,fs);
    sonoAux = sonogram(lfpAux,fs,0.9,2.5);
    
    N = length(lfpAux);
    stmSub = [];
    for cc = 1:numel(Conditions)
        trgSub = Conditions{cc}.Triggers;
        disp(Conditions{cc}.name)
        if isrow(trgSub)
            trgSub = trgSub';
        end
        
        stmSub =[stmSub; trgSub, trgSub + round(5.35*fs)];
        if strcmp(Conditions{cc}.name,'MechControl')
            mecSub = [trgSub, trgSub+ round(5.35*fs)];
        end
    end
    stmSub = round(stmSub);
    mecSub = round(mecSub);
    stimTrace = false(1,N);
    mechTrace = false(1,N);
    for ct = 1:size(stmSub)
        stimTrace(stmSub(ct,1):stmSub(ct,2)) = true;
    end
    for cmt = 1:size(mecSub)
        mechTrace(mecSub(cmt,1):mecSub(cmt,2)) = true;
    end
    silIdx = ~stimTrace;
    %% Computing the spectrum from the desired segments.
    sponSpect = getMeanSpectrumFromSegments(silIdx,sonoAux,fs);
    mechSpect = getMeanSpectrumFromSegments(mechTrace,sonoAux,fs);
    
    SAL_MECH{cf} = {[sonoAux.FrequencyAxis;mechSpect]};
    SAL_SPCT(cf,1) = {[sonoAux.FrequencyAxis;abs(sponSpect)]};
    SAL_SPCT(cf,2) = {[sonoAux.FrequencyAxis;angle(sponSpect)]};
    
    subplot(1,2,1);plot(sonoAux.FrequencyAxis,...
        20*log10(abs(sponSpect)/length(sponSpect)),'DisplayName',dirFiles(cf).name)
    if cf == 1
        hold on
    end
    subplot(1,2,2);plot(sonoAux.FrequencyAxis,unwrap(angle(sponSpect)))
    if cf == 1
        hold on
    end
    SAL_LFP(cf,1) = {EEG.data};
    SAL_LFP(cf,2) = {EEG.header.SamplingFrequency};
    SAL_LFP(cf,3) = {dirFiles(cf).name};
    SAL_LFP(cf,4) = {Conditions};
end
subplot(1,2,1);title('Power spectrum for Saline LFPs');xlabel('Frequency [Hz]')
ylabel('dB')
subplot(1,2,2);title('Angular spectrum for Saline LFPs');xlabel('Frequency [Hz]')
ylabel('Degrees/Radians')
