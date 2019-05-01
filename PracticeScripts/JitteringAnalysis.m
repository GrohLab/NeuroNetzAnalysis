%% Constant values
% Prefixes
k = 1e3; % Kilo
m = 1e-3; % milli 

%% Select data file
defaultPath = 'E:\Data\Jittering'; % Path where you would like to open a windows explorer window to select your data file
[fileName, filePath] =...
    uigetfile({'*.smr';'*.smrx';'*.mat'},'Select a data file',defaultPath);
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
    case '.mat'
        fprintf(1,'\nMatlab file selected.\n')
    otherwise
        fprintf(1,'\nFile type not recognized. Please provide a valide file\n')
        return
end
%% Loading the variables to extract all triggers
% Loading the triggers (light/laser and piezo pulses)
varsInFile = load(fullfile(filePath,fileName));
fNames = fieldnames(varsInFile);
chanIdx = startsWith(fNames,'chan');
headIdx = startsWith(fNames,'head');
headSub = find(headIdx);
fsArray = zeros(1,nnz(headIdx));
for chead = 1:nnz(headIdx)
    fsArray(chead) = 1e6/(varsInFile.(fNames{headSub(chead)}).sampleinterval);
    newField = erase(varsInFile.(fNames{headSub(chead)}).title,' ');
    data.(newField) = varsInFile.(strrep(fNames{headSub(chead)},'head','chan'));
end
%%
timeLapse = [50*m,150*m]; % Time window surrounding the trigger [time before, time after] in seconds
fs = 20*k; % Sampling frequency normally 20 kHz

% Creation of the discrete and the continuous stacks:
%                                        First spikeing times  trigger  timeW     FS  FSt  N spiking times          continuous data. 
[discreteStack, contStack] = getStacks(spikeFindingData.spikes,upT,'on',timeLapse,fs,fs,spikeFindingData2.spikes,EEG.data);
%                                                               indx                   cell of indexes | INDX,  signal1, signal2,...
[Ne, Nt, Na] = size(discreteStack);
kIdx = false(1,Na);
binSz = 10; % milliseconds or seconds
ERASE_kIDX = false;
[PSTH, trig, sweeps, tx_psth] = getPSTH(discreteStack, timeLapse, kIdx, binSz, fs);
[triggeredAverageSignals, signalVariation, tx_tav] = getTriggeredAverage(contStack, kIdx, timeLapse);
[relativeSpikeTimes, tx_raster] = getRasterFromStack(discreteStack, kIdx, timeLapse, fs, ERASE_kIDX);

figure;plot(tx_psth,PSTH(1,:),tx_psth,PSTH(2,:))
plotRaster(relativeSpikeTimes, timeLapse, fs, '130626_p3');
plotRasterAsText(relativeSpikeTimes, timeLapse, fs, '130626_p3')