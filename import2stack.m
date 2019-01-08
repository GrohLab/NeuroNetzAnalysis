function iOk = import2stack(EphysPath)
%IMPORT2STACK uses the file system from Anton to create a single
%*analysis.mat file.
%   The function accepts the path to the file system and creates the
%   non-existent *analysis.mat files. 
if ~exist(EphysPath,'dir')
    fprintf('The given path does not exist.\n')
    iOk = -1;
else
    inDirStruct = dir(EphysPath);
    dirContents = extractfield(inDirStruct,'name');
    searchDirs = {'continuousEphys', 'LFP', 'LightDetail', 'Whisker',...
        'DiscreteData.mat','EmptyHead.mat'};
    [~, idxDir] = intersect(dirContents,searchDirs);
    if length(idxDir) == 6
        if ~isfolder(fullfile(EphysPath,'AnalysisMatFiles'))
            mkdir(EphysPath,'AnalysisMatFiles')
        end
        load(fullfile(EphysPath,'DiscreteData.mat'),'DiscreteData')
        load(fullfile(EphysPath,'EmptyHead.mat'),'ExampleHead')
        load(fullfile(EphysPath,'..','ephys_database.mat'),'RecDB','ExpDB')
        Nfiles = length(DiscreteData);
        for cf = 1:Nfiles
            expDD = DiscreteData(cf);
            expName = RecDB.Properties.RowNames{cf};
            fs = 2e4;fsLFP = 1e3;
            anaFile = fullfile(EphysPath,'AnalysisMatFiles',expName,...
                    [expName,'analysis.mat']);
            if ~exist(anaFile,'file')
                % Create the analysis file and directory if it doesn't
                % exist
                if ~exist(getParentDir(anaFile,1),'dir')
                    mkdir(getParentDir(anaFile,1))
                end
                %% Continuous signals
                LFPprobeDepth = ExpDB{{RecDB.AnimalName{cf}}, 'LfpCoord'}(3); %#ok<IDISVAR,USENS>
                [LFP, whisker] = loadLFPAndWhisker(LFPprobeDepth,expName,EphysPath);
                load(fullfile(EphysPath,searchDirs{1},[expName,'.mat']),'fR','rR')
                lightSignal = geEphystStimPeriods(expDD,fs,fs,'l');
                puffSignal = getStimPeriods(expDD,fs,fs,'p');
                touchSignal = getStimPeriods(expDD,fs,fs,'t');
                whiskingSignal = getStimPeriods(expDD,fs,fs,'w');
                poleSignal = getStimPeriods(expDD,fs,fs,'pl');
                groomingSignal = getStimPeriods(expDD,fs,fs,'g');
                excludeSignal = getStimPeriods(expDD,fs,fs,'x');
                %% Spikes information 'spikesFindingData'
                [~,sp,~] = getStimPeriods(expDD,fs,fs,'s');
                minISI = 1000*(min(diff(sp))/fs); % Minimal ISI in ms.
                ppms = fs/1000;
                %% Building the analysis structures
                spikeFindingData = struct('thresh',[],'minISI',minISI,...
                    'spikes',sp,'ppms',ppms,'timestamp',[]); %#ok<*NASGU>
                filteredResponse = createStructure(fR,fs,ExampleHead);
                RawResponse = createStructure(rR,fs,ExampleHead);
                EEG = createStructure(LFP,fsLFP,ExampleHead);
                Triggers = struct('whisker',whiskingSignal,'light',lightSignal,...
                    'puff',puffSignal,'touch',touchSignal,...
                    'whisking',whisker,'pole',poleSignal,...
                    'grooming',groomingSignal,'exclude',excludeSignal);
                notes = RecDB.PhysioNucleus(cf);
                Conditions = {};
                save(analysisFile,...
                    'notes','RawResponse','Triggers','filteredResponse','EEG',...
                    'spikeFindingData','Conditions')
            else
                fprintf('%s Analysis file exists. Skipping experiment...\n',...
                    expName)
            end
            iOk = 1;
        end
    else
        fprintf('Some folders/files are missing!\n')
        iOk = -1;
    end
end

end

function [stimPeriods,stimStart,stimRF] = getStimPeriods(dd,fs,fs2,stimString)
% 'Down-sampling' the time indeces for the stimulus/triggers
fact = fs2/fs;
stimStart = [];
stimLength = [];
switch stimString
    case 's'
        stimStart = round(dd.Spikes*fact);
        stimLength = zeros(size(stimStart));
    case 'w'
        stimStart = round(dd.WhiskingStart*fact);
        stimLength = round(dd.WhiskingLength*fact);
    case 'l'
        stimStart = round(dd.LightStart*fact);
        stimLength = round(dd.LightLength*fact);
    case 'p'
        stimStart = round(dd.PuffStart*fact);
        stimLength = round(dd.PuffLength*fact);
    case 't'
        stimStart = round(dd.TouchStart*fact);
        stimLength = round(dd.TouchLength*fact);
    case 'pl'
        if ~isempty(dd.Sections.Pole)
            stimStart = round(dd.Sections.Pole(:,1)*fact);
            stimLength = round((dd.Sections.Pole(:,2)-dd.Sections.Pole(:,1))...
                * fact);
        end
    case 'g'
        if ~isempty(dd.Sections.Grooming)
            stimStart = round(dd.Sections.Grooming(:,1)*fact);
            stimLength = round((dd.Sections.Grooming(:,2)-dd.Sections.Grooming(:,1))...
                * fact);
        end
    case 'x'
        if ~isempty(dd.Sections.Exclude)
            stimStart = round(dd.Sections.Exclude(:,1)*fact);
            stimLength = round((dd.Sections.Exclude(:,2)-dd.Sections.Exclude(:,1))...
                * fact);
        end
    otherwise
        stimPeriods = dd.LengthTime;
        stimStart = 0;
        fprintf('No stimuli recognized...')
        return;
end
stimPeriods = false(1,round(dd.LengthTime*fs2));
for counter = 1:length(stimStart)
    stimPeriods(...
        stimStart(counter):stimStart(counter)+stimLength(counter)) = true;
end
stimRF = [stimStart',stimLength'];
end

function [LFP,whisker] = loadLFPAndWhisker(LFPprobeDepth,ExpName,EphysPath)
% LFP is recorded in 16 linear channels along cortex. The order from white
% matter to pia is:
LFPsort=[6, 11, 3, 14, 1, 16, 2, 15, 5, 12, 4, 13, 7, 10, 8, 9];
% The individual channels depth is:
LFPdepth=LFPprobeDepth - (50:100:1550);
% Consequently for ~L5 lfp at 850 ?m depth, when LFPprobeDepth
% is 1600 would be:
LFPchIdx = LFPsort(find(LFPdepth < 900, 1));
% Loading only the L5 LFP from the 16 channels file
LFP_L5 = sprintf('LFPch%d',LFPchIdx);
try
    LFP = load([EphysPath, 'LFP\', ExpName, '.mat'],LFP_L5);
    LFP = struct2array(LFP);
catch
    LFP = [];
end
try
    whisker = load([EphysPath, 'Whisker\', ExpName, '.mat'],...
        'WhiskerAngle');
    whisker = struct2array(whisker);
catch
    whisker = [];
end
end

function varStruct = createStructure(chan,fs,exHead)
if ~isempty(chan)
    auxHead = getChannelHeader(exHead,fs,chan);
    varStruct = struct('data',chan,'header',auxHead);
else
    varStruct = struct('data',[],'header',exHead);
end
end

function head = getChannelHeader(exHead,fs,chan)
head = exHead;
head.min = min(chan);
head.max = max(chan);
head.SamplingFrequency = fs;
head.npoints = numel(chan);
end