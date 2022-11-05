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
directory = 'E:\LFP_Frequency_Analysis\Juxta\';
inDirs = dir(directory);
ccon = 1;
ovrlp = 0.95; % percentage
winSz = 5; % seconds
spectra = {};
ITPC = {};
while ccon <= numel(inDirs)
    if inDirs(ccon).name == "." || inDirs(ccon).name == ".." || ...
            ~inDirs(ccon).isdir
        inDirs(ccon) = [];
        disp('File or dot(s) found')
        continue;
    end
    fprintf('Processing %s folder...\n',inDirs(ccon).name)
    conNames(ccon) = {inDirs(ccon).name};
    lfpStack = [];
    dirFiles = dir(fullfile(inDirs(ccon).folder,inDirs(ccon).name,'*analysis.mat'));
    
    for cf = 1:length(dirFiles)
        %% Loading and a small preprocessing
        load(fullfile(dirFiles(cf).folder,dirFiles(cf).name),...
            'EEG','Conditions')
        fs = EEG.header.SamplingFrequency;
        desFs = 1e3;
        lfpAux = EEG.data;
        if ~isa(lfpAux,'double') && iscolumn(lfpAux)
            lfpAux = double(lfpAux)';
        end
        %% Resampling the original signal
        N_orig = length(lfpAux);
        tx = (0:N_orig-1)*(1/fs);
        desTx = 0:1/desFs:tx(end);
        lfpTS = timeseries(lfpAux,tx);
        lfpTS = lfpTS.resample(desTx);
        lfpAux = squeeze(lfpTS.Data)';
        lfpAux = brainwaves(lfpAux,desFs,{'alpha',1e-2,200});
        lfpAux = iir50NotchFilter(lfpAux,desFs);
        N_subsamp = length(lfpAux);
        
        %% Transforming the triggers' time points
        stmSub = [];
        for cc = 1:numel(Conditions)
            trgSub = Conditions{cc}.Triggers;
            trgSub = round(trgSub .* (fs\desFs));
            
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
            getStacks(false(1,N_subsamp),...
            mecSub,'on',[1,6],desFs,desFs,[],lfpAux);
        lfpStack = cat(3,lfpStack,zeros(1,size(stackAux,2),1),stackAux);
        %% Trying the sonogram once again
        sonoCorr = (winSz*desFs)/2;
        stmIdx = StepWaveform.subs2idx(stmSub - sonoCorr,N_subsamp);
        spnIdx = ~stmIdx;
        mecIdx = StepWaveform.subs2idx(mecSub - sonoCorr,N_subsamp);
        sonoStruct = sonogram(lfpAux,desFs,ovrlp,winSz);
        [spnSpect,spnITPC] = getMeanSpectrumFromSegments(spnIdx, sonoStruct, desFs);
        [mecSpect,mecITPC] = getMeanSpectrumFromSegments(mecIdx, sonoStruct, desFs);
        spectra(ccon,1:2) = {spnSpect,mecSpect};
        ITPC(ccon,1:2) = {spnITPC,mecITPC};
        
    end
    %% Save results
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
    
    N_orig = length(lfpAux);
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
    stimTrace = false(1,N_orig);
    mechTrace = false(1,N_orig);
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
