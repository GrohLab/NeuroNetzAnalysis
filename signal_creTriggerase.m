function [dStruct, cStruct] =...
    signal_creTriggerase(configStruct,analysisFilePath)
%SIGNAL_TRIGGERASE cuts the signals in the analysis file specified by
%_analysisFileName_ given the trigger signal in the configuration structure
% (_configStruct_). 
%    This function does not need the exclusion nor the neglection
%    information contained in the configuration structure; it creates the
%    stacks regardless of this information.

% Helpful function handlers 
expRng = @(x) x(1):x(2);

%% Validation of the structure
if isempty(configStruct) || isempty(configStruct.Trigger) ||...
        isempty(configStruct.Trigger.Name)
    fprintf('The stack cannot be built without a trigger.\n')
    fprintf('Please check this issue and try again creating a valid ')
    fprintf('configuration structure.\n')
    dStruct = [];
    cStruct = [];
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
    twfObj = StepWaveform(Sig.(sigIDs{tnIdx}),fs);
    tSubs = twfObj.Triggers;
    dsIdx = structfun(@islogical,Sig);
    Nds = sum(dsIdx);
    lenFlds = structfun(@length,Sig);
    discEvents = false(Nds,ceil(mean(lenFlds(dsIdx))));
    cds = 1;
    for cf = 1:numel(sigIDs)
        if islogical(Sig.(sigIDs{cf}))
            discEvents(cds,:) = Sig.(sigIDs{cf});
            cds = cds + 1;
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
ndName = sigIDs{~dsIdx};
sigIDs(~dsIdx) = [];
sigIDs = cat(isrow(sigIDs)*2 + iscolumn(sigIDs)*1,{'spikes'},sigIDs);
discEvents = cat(1,StepWaveform.subs2idx(spSubs,size(discEvents,2)),...
    discEvents);
discEvents = [discEvents(strcmpi(sigIDs,configStruct.Trigger.Name),:);...
    discEvents];
discEvents([false;strcmpi(sigIDs,configStruct.Trigger.Name)],:) = [];
sigIDs = cat(isrow(sigIDs)*2 + iscolumn(sigIDs)*1,...
    {configStruct.Trigger.Name},sigIDs);
sigIDs(...
    cat(isrow(sigIDs)*2 + iscolumn(sigIDs)*1,...
    false,...
    strcmpi(sigIDs(2:end),configStruct.Trigger.Name))) = [];


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
if numel(Sig) == 1
    WHI = Sig.(ndName);
else
    fprintf('It is unclear where to take the whisker signal from.\n')
    WHI = zeros(size(LFP.data),'single');
end
contEvents = cat(1,LFP.data,WHI);
conIDs = {'LFP','Whisking'};
NsC = length(contEvents);

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
tSubsC = round(tSubs .* (fs(2)/fs(1)));
continuouStack = zeros(NeC,NtC,Na);

%% Construct stacks: cre_Triggerase
for cap = 1:size(tSubs,1)
    onf = configStruct.Trigger.Edge*1 + ~configStruct.Trigger.Edge*2;
    ctSubs = tSubs(cap,onf);
    ctSubsC = tSubsC(cap,onf);
    segSubs = [ctSubs - prevSamples, ctSubs + postSamples];
    normSubs = segSubs;
    if segSubs(1) < 1
        printf('The viewing window is out of the signal range.\n')
        fprintf('%d samples required before the signal. Omitting...\n',...
            -segSubs(1))
        segSubs(1) = 1;
        normSubs = Nt - segSubs(2) + 1:Nt;
    elseif segSubs(2) > Ns
        fprintf('The viewing window is out of the signal range.\n')
        fprintf('%d samples required after the signal. Omitting...\n',...
            segSubs(2))
        segSubs(2) = Ns;
        normSubs = 1:diff(segSubs)+1;
    end
    if diff(fs) %% Different signal lengths
        segSubsC = [ctSubsC - prevSamplesC, ctSubsC + postSamplesC];
        normSubsC = segSubsC;
        if segSubsC(1) < 1
            printf('The viewing window is out of the signal range.\n')
            fprintf('%d samples required before the signal. Omitting...\n',...
                -segSubsC(1))
            segSubsC(1) = 1;
            normSubsC = NtC - segSubsC(2) + 1:NtC;
        elseif segSubsC(2) > NsC
            fprintf('The viewing window is out of the signal range.\n')
            fprintf('%d samples required after the signal. Omitting...\n',...
                segSubsC(2))
            segSubsC(2) = NsC;
            normSubsC = 1:diff(segSubsC)+1;
        end
    else
        segSubsC = segSubs;
        normSubsC = normSubs;
    end
end
discreteStack(:,expRng(normSubs),cap) = discEvents(:,expRng(segSubs));
continuouStack(:,expRng(normSubsC),cap) = contEvents(:,expRng(segSubsC));
%% Output structures:
dStruct = struct('Stack',discreteStack,'SignalIDs',sigIDs);
cStruct = struct('Stack',continuouStack,'SignalIDs',conIDs);
end