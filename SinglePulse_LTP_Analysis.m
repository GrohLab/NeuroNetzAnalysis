%% LTP
% Loading the necessary files (spike times, 
dataDir = 'E:\Data\VPM\LTP\190701_LTP_3700_1500_1520';
figureDir = fullfile(dataDir,'Figures\');
if ~loadTriggerData(dataDir)
    fprintf(1,'Not possible to load all the necessary variables')
    return
end
%% Constructing the helper 'global' variables
% Number of conditions
Nccond = numel(Conditions);
% Number of total samples
Ns = min(structfun(@numel,Triggers));
% Total duration of the recording
Nt = Ns/fs;
% Useless clusters (labeled as noise or they have very low firing rate)
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
bads = find(badsIdx);
totSpkCount = cellfun(@numel,sortedData(:,2));
clusterSpikeRate = totSpkCount/Nt;
silentUnits = clusterSpikeRate < 0.1;
bads = union(bads,find(silentUnits));
goods = setdiff(1:size(sortedData,1),bads);
badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));
% Logical spike trace for the first good cluster
spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
% Subscript column vectors for the rest good clusters
spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods(2:end),2),...
    'UniformOutput',false);
% Number of good clusters 
Ncl = numel(goods);
% Redefining the stimulus signals from the low amplitude to logical values
whStim = {'piezo','whisker'};
cxStim = {'laser','light'};
lfpRec = {'lfp','s1','cortex','s1lfp'};
trigNames = fieldnames(Triggers);
numTrigNames = numel(trigNames);
ctn = 1;
while ctn <= numTrigNames 
    if contains(trigNames{ctn},whStim,'IgnoreCase',true)
        whisker = Triggers.(trigNames{ctn})(1:Ns);
    end
    if contains(trigNames{ctn},cxStim,'IgnoreCase',true)
        laser = Triggers.(trigNames{ctn})(1:Ns);
    end
    if contains(trigNames{ctn},lfpRec,'IgnoreCase',true)
        LFP = Triggers.(trigNames{ctn})(1:Ns);
    end
    ctn = ctn + 1;
end
mObj = StepWaveform(whisker,fs);
mSubs = mObj.subTriggers;
piezo = mObj.subs2idx(mSubs,Ns);
lObj = StepWaveform(laser,fs);
lSubs = lObj.subTriggers;
laser = lObj.subs2idx(lSubs,Ns);
mObj.delete;lObj.delete;
continuousSignals = {piezo;laser;LFP};
%% User prompt for relevant information:
% Time lapse, bin size, and spontaneous and response windows
promptStrings = {'Viewing window (time lapse) [s]:','Response window [s]',...
    'Bin size [s]:'};
defInputs = {'0.1, 0.1', '0.002, 0.05', '0.01'};
answ = inputdlg(promptStrings,'Inputs', [1, 30],defInputs);
if isempty(answ)
    fprintf(1,'Cancelling...\n')
    return
else
    timeLapse = str2num(answ{1}); %#ok<*ST2NM>
    if numel(timeLapse) ~= 2
        timeLapse = str2num(inputdlg('Please provide the time window [s]:',...
            'Time window',[1, 30], '0.1, 0.1'));
        if isnan(timeLapse) || isempty(timeLapse)
            fprintf(1,'Cancelling...')
            return
        end
    end
    responseWindow = str2num(answ{2});
    binSz = str2double(answ(3));
end
fprintf(1,'Time window: %.2f - %.2f ms\n',timeLapse(1)*1e3, timeLapse(2)*1e3)
fprintf(1,'Response window: %.2f - %.2f ms\n',responseWindow(1)*1e3, responseWindow(2)*1e3)
fprintf(1,'Bin size: %.3f ms\n', binSz*1e3)
spontaneousWindow = -flip(responseWindow);
%% Condition triggered stacks
condNames = arrayfun(@(x) x.name,Conditions,'UniformOutput',false);
% Choose the conditions to create the stack upon
[chCond, iOk] = listdlg('ListString',condNames,'SelectionMode','single',...
    'PromptString','Choose the condition to look at:');
if ~iOk
    fprintf(1,'Cancelling...\n')
    return
end

% Select the onset or the offset of a trigger
fprintf(1,'Condition ''%s''\n', Conditions(chCond).name)
onOffStr = questdlg('Trigger on the onset or on the offset?','Onset/Offset',...
    'on','off','Cancel','Onset');
if strcmpi(onOffStr,'Cancel')
    fprintf(1,'Cancelling...\n')
    return
end

% Constructing the stack out of the user's choice
[dst, cst] = getStacks(spkLog,Conditions(chCond).Triggers,onOffStr,...
    timeLapse,fs,fs,spkSubs,continuousSignals);

[Ne, Nt, NTa] = size(dst);
% Computing the time axis for the stack
tx = (0:Nt)/fs - timeLapse(1);
%% LTP Exclusive: Control and post-induction conditions
% Boolean flags indicating which trigger belongs to which condition (delay
% flags)
% Determining if the gap between stimuli are big enough to separate the
% conditions
avGap = mean(diff(Conditions(chCond).Triggers(:,1),1,1));
[gapVal, timeGapSub] = max(diff(Conditions(chCond).Triggers(:,1),1,1));
errGap = log(avGap) - log(gapVal);
if abs(errGap) > 1
    % If the gap is bigger in one order of magnitude, we consider it to be
    % another condition (Control / Post-induction)
    condFlags = false(NTa, 2);
else
    % Otherwise, it is only one condition
    condFlags = true(NTa, 1);
end
Na = sum(condFlags,1);
%% Counting spikes in given windows
sponActStackIdx = tx >= spontaneousWindow(1) & tx <= spontaneousWindow(2);
respActStackIdx = tx >= responseWindow(1) & tx <= responseWindow(2);

