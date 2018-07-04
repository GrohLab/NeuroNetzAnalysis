%% Create the .mat file from the .smrx for Sailaja's Data
% "C:\Users\neuro\Documents\MATLAB\Data\M104_C5_Mech_5mW.smrx"
% Run this section or use the DataAnalysis GUI to create both the .mat and
% the *analysis.mat files.

%% Regardless of the DataAnalysis GUI, you should give the information in this section
baseName = 'M104_C5_Mech_5mW';
smrxFileName = [baseName, '.smrx'];
filePath = 'C:\Users\neuro\Documents\MATLAB\Data\';
%% If you used DataAnalysis GUI, you shouldn't run this section
importSMR(smrxFileName, filePath, 256);
%% Use UMS2k for the spike extraction. You should run this setion once
matFileName = fullfile(filePath,[baseName,'.mat']);
load(matFileName, 'chan6') % This channel is the pressure signal. IMPORTANT
readingChannels = 9; % This is the filtered response channel. IMPORTANT
experimentObject =...
    UMSDataLoader(matFileName, readingChannels);
disp(experimentObject)
%%
experimentObject.UMS2kPipeline;
experimentObject.getSpikeTimes;
experimentObject.plot;
%% Save the spikes if you agree
experimentObject.saveSpikeTimes;
%% Triggered stack. You can run this section more than once
sortedFile = fullfile(filePath,[baseName,'sorted.mat']);
% sortedFile = fullfile(filePath,[baseName,'UMS_analysis.mat']);
if exist(sortedFile,'file') && experimentObject.Ns > 100
    experimentObject.loadSpikeTimes;
elseif exist(sortedFile,'file')
    experimentObject.loadSpikeTimes(sortedFile);
else
    fprintf('The sorted file doesn''t exist... yet\n')
end

%% Initialization for the loop
analysisFileName =...
    fullfile(filePath,[baseName,'analysis.mat']);
load(analysisFileName,'Conditions','EEG','Triggers')
Ncon = numel(Conditions);
expStack = cell(1,Ncon);
LFPstack = expStack;
MechStack = LFPstack;
timeLapse = [2, 7];       % Time before the trigger, time after the trigger
cellType = 'other';
binningTime = 100e-3;    % 10 
fs = experimentObject.SamplingFrequency;

%% Looping throught the found conditions.
for ccon = 1:Ncon
    % Create stacks
    [auxExp, auxCont] = getStacks(...      CREATE STACK
        experimentObject.SpikeTimes,...     Spike times
        Conditions{ccon}.Triggers, 'on',...    Onset of triggers
        timeLapse,...                       
        fs,fs,[],...   No extra events
        EEG.data,chan6,Triggers.whisker,Triggers.light);
    expStack(ccon) = {auxExp};
    LFPstack(ccon) = {squeeze(auxCont(1,:,:))};
    MechStack(ccon) = {squeeze(auxCont(2,:,:))};
    auxMech = squeeze(auxCont(3,:,:));
    auxLight = squeeze(auxCont(4,:,:));
    fprintf('Finished creating the stack\nCreating plot...\n')
    [~,FigIDNW] =...
        plotTriggeredEvents(...             PLOT TRIGGERS
        expStack{ccon},...                  Current condition
        LFPstack{ccon},...                  LFP
        MechStack{ccon},...                 Mechanical pressure
        timeLapse,...                       Time before, time after
        cellType,...                        'POm', 'VPM', or 'other'
        binningTime,...                     Bin size in seconds
        fs,fs);
    fprintf('Plot created!\n')
    title(Conditions{ccon}.name)
    hold on
    %% Adding the continuous triggered average.
    % Pressure signal
    fprintf('Adding the TTL mechanical stimulus\n')
    [pressureSign,~] =...
        getTriggeredAverage(MechStack{ccon},false(1,size(MechStack{ccon},2)));
    % Mechanical action TTL
    [meanMech, ~] = getTriggeredAverage(auxMech,false(1,size(auxMech,2)));
    tx = 0:1/fs:(length(meanMech)-1)/fs;
    tx = tx - timeLapse(1);
    meanMech = equalizeAmplitudeTo(pressureSign,meanMech);
    plot(tx,meanMech,'LineWidth',0.6,'Color',[37, 154, 3]/255)
    text(tx(end),mean(double(meanMech)),'Pressure','FontWeight','bold',...
        'HorizontalAlignment','right','Color',[37, 154, 3]/255)
    fprintf('Adding the TTL laser stimulus\n')
    
    % Laser TTL
    [meanLight, ~] = getTriggeredAverage(auxLight,false(1,size(auxLight,2)));
    meanLight = equalizeAmplitudeTo(pressureSign, meanLight);
    plot(tx,meanLight,'LineWidth',0.6,'Color',[0, 138, 230]/255)
    text(tx(end),mean(double(meanLight)),'Light','FontWeight','bold',...
        'HorizontalAlignment','right','Color',[0, 138, 230]/255)
    % Add the printing fixed options with the file name without latex
    % interpretation
end

