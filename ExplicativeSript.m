%% Create the .mat file from the .smrx for Sailaja's Data
% "C:\Users\neuro\Documents\MATLAB\Data\M104_C5_Mech_5mW.smrx"
% Run this section or use the DataAnalysis GUI to create both the .mat and
% the *analysis.mat files.

baseName = 'M104_C5_Mech_5mW';
smrxFileName = [baseName, '.smrx'];
filePath = 'C:\Users\neuro\Documents\MATLAB\Data\';
importSMR(smrxFileName, filePath, 256);
%% Use UMS2k for the spike extraction
matFileName = fullfile(filePath,[baseName,'.mat']);
load(matFileName, 'chan6') % This line is not nec
readingChannels = 9; %%%%%%%%%%%%%%%%%%%%%%% Fileted response channel! 
experimentObject =...
    UMSDataLoader(matFileName, readingChannels);
disp(experimentObject)
experimentObject.UMS2kPipeline;
experimentObject.getSpikeTimes;
experimentObject.plot;
experimentObject.saveSpikeTimes;
%% Triggered stack
analysisFileName =...
    fullfile(filePath,[baseName,'analysis.mat']);
load(analysisFileName,'Conditions','EEG','Triggers')
Ncon = numel(Conditions);
expStack = cell(1,Ncon);
LFPstack = expStack;
MechStack = LFPstack;
timeLapse = [1, 6];       % Time before the trigger, time after the trigger
cellType = 'other';
binningTime = 100e-3;    % 10 
fs = experimentObject.SamplingFrequency;

for ccon = 1:Ncon
    % Create stack
    [auxExp, auxLFP, auxMech] = getStack(...      CREATE STACK
        experimentObject.SpikeTimes,...     Spike times
        Conditions{ccon}.Triggers, 'on',...    Onset of triggers
        timeLapse,...                       
        fs,...
        EEG.data,chan6,[],fs);                          % No extra events
    expStack(ccon) = {auxExp};
    LFPstack(ccon) = {auxLFP};
    MechStack(ccon) = {auxMech};
    [~,FigIDNW] =...
        plotTriggeredEvents(...             PLOT TRIGGERS
        expStack{ccon},...                  Current condition
        LFPstack{ccon},...                  LFP
        MechStack{ccon},...                 Mechanical pressure
        timeLapse,...                       Time before, time after
        cellType,...                        'POm', 'VPM', or 'other'
        binningTime,...                     Bin size in seconds
        fs,fs);
    title(Conditions{ccon}.name)
    clearvars aux*
    [~,auxMech,auxLight] = getStack(...
        experimentObject.SpikeTimes,...
        Conditions
end