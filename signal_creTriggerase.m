function [discreteStack,continuouStack] =...
    signal_creTriggerase(configStruct,analysisFilePath)
%SIGNAL_TRIGGERASE cuts the signals in the analysis file specified by
%_analysisFileName_ given the trigger signal in the configuration structure
% (_configStruct_). 
%    This function does not need the exclusion nor the neglection
%    information contained in the configuration structure; it creates the
%    stacks regardless of this information.
%% Validation of the structure
if isempty(configStruct) || isempty(configStruct.Trigger) ||...
        isempty(configStruct.Trigger.Name)
    fprintf('The stack cannot be built without a trigger.\n')
    fprintf('Please check this issue and try again creating a valid ')
    fprintf('configuration structure.\n')
    discreteStack = [];
    continuouStack = [];
    return
end
%% Loading and acconditioning data from the analysis file. NOT GENERALIZED
% Reading the sampling frequency of the experiment
fInfo = load(analysisFilePath,'filteredResponse');
fs = fInfo.filteredResponse.header.SamplingFrequency;
Ns = fInfo.filteredResponse.header.npoints;
% Loading the spiking data from UMS. If the UMS pipeline was not ran, the
% spiking data is from the threshold method. BEWARE of the UMSDataLoader
% instance created in this part of the code!
UMS = load(analysisFilePath,'UMSSpikeStruct');
if ~isempty(UMS)
    UMS = UMS.UMSSpikeStruct;
    fs = UMS.params.Fs;
    spikeLoader = UMSDataLoader();
    spikeLoader.changeUMSStructure(UMS)
    spikeLoader.getSpikeTimes
    spSubs = spikeLoader.SpikeTimes;
else
    UMS = load(analysisFilePath,'spikeFindingData');
    UMS = UMS.spikeFindingData;
    spSubs = UMS.spikes;
end

% Loading events from the experiment. If it is not existing, then the
% algorithm needs to load ConditionsTest as a last attempt to continue.
% BEWARE of the StepWaveform class static method used in this part of the
% code!
Sig = load(analysisFilePath,'Triggers');
if ~isempty(Sig)
    Sig = Sig.Triggers;
    sigIDs = fieldnames(Sig);
    tnIdx = strcmpi(sigIDs,configStruct.Trigger.Name);
    twfObj = StepWaveform(Sig.(sigIDs(tnIdx)),fs);
    tSubs = twfObj.Triggers;
    dsIdx = structfun(@islogical,Sig);
    Nds = sum(dsIdx);
    lenFlds = structfun(@length,Sig);
    discEvents = false(Nds,mean(lenFlds(dsIdx)));
    for cf = 1:numel(sigIDs)
        if islogical(Sig.(sigIDs(cf)))
            discEvents(cf) = Sig.(sigIDs(cf));
        end
    end
else
    Sig = load(analysisFilePath,'ConditionsTest');
    Sig = Sig.ConditionsTest;
    Sig = cell2mat(Sig);
    getname = @(x) x.name;
    sigIDs = arrayfun(getname,CT2,'UniformOutput',false);
    tnIdx = strcmpi(sigIDs,configStruct.Trigger.Name);
    tSubs = Sig(tnIdx).Triggers;
    discEvents = false(numel(sigIDs),Ns);
    for cf = 1:numel(sigIDs)
        discEvents(cf,:) = StepWaveform.subs2idx(Sig(cf).Triggers,Ns);
    end
end
sigIDs(~dsIdx) = [];
sigIDs = cat(isrow(sigIDs)*2 + iscolumn(sigIDs)*1,{'spikes'},sigIDs);
discEvents = cat(1,StepWaveform.subs2idx(spSubs,size(discEvents,2)),...
    discEvents);
% Loading the local field potential signal (sampled at 1 kHz for Toni's
% data). There is no alternative if this variable is missing.
LFP = load(analysisFilePath,'EEG');
LFP = LFP.EEG;
if fs ~= LFP.header.SamplingFrequency
    fs(2) = LFP.header.SamplingFrequency;
else
    fs = repmat(fs,2,1);
end
% Get the whisking signal
wIdx = strcmpi(sigIDs,'whisking');
if numel(Sig) > 1
    WHI = Sig.(sigIDs(wIdx));
else
    fprintf('It is unclear where to take the whisker signal from.\n')
    WHI = zeros(size(LFP.data),'single');
end
contEvents = cat(1,LFP.data,WHI);
conIDs = {'LFP','Whisking'};

%% Construct stacks: Initialization
%   Preallocation of the discrete stack:
timeSpan = configStruct.ViewWindow;
prevSamples = ceil(timeSpan(1) * fs(1));
postSamples = ceil(timeSpan(2) * fs(1));
% Total number of samples per trial
Nt = prevSamples + postSamples + 1;
% Total number of discrete conditioning variables
Ne = size(discEvents,1);
% Total number of trigger points
Na = size(tSubs,1);
discreteStack = false(Ne,Nt,Na);

%   Preallocation of the continuous stack:
prevSamplesC = ceil(timeSpan(1) * fs(2));
postSamplesC = ceil(timeSpan(2) * fs(2));
% Number of time points
NtC = prevSamplesC + postSamplesC + 1;
% Total number of continuous conditioning variables.
NeC = size(contEvents,1);
continuouStack = zeros(NeC,NtC,Na);

%% Construct stacks: cre_Triggerase

for cap = 1:size(tSubs,1)
    
end


end