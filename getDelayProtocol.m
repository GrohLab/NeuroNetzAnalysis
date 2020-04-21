function [Conditions, Triggers] = getDelayProtocol(expFolder)
% GETDELAYPROTOCOL searches specifically for the stimulation conditions for
% the jittering project. The delays can vary and the whisker stimulation
% should always have a control pulse, as well as the cortex excitation
% using laser pulses. WARNING! No frequency stimulation should be given to
% this function.

%% Gathering the information about the experiment

checkSignal = @(x,y) contains(x,y,'IgnoreCase',true);
findUnpairedPulse = @(x) cat(1,false,reshape(diff(sort(x)) == 0,numel(x)-1,1));
Conditions = struct('name',{},'Triggers',{});
Triggers = struct();
if ~loadTriggerData(expFolder)
    fprintf(1,'No condition extracted.\n')
    return
end
%{
smrxFile = dir(fullfile(expFolder,'*.smrx'));
[~,expName,~] = fileparts(smrxFile(1).name);
if isempty(smrxFile)
    warning('There is no experiment in this folder...\n')
    return
end
if exist(fullfile(expFolder,[expName,'analysis.mat']),'file')
    answ = questdlg(...
        'The analysis file for this experiment exists already. Would you like to overwrite it?',...
        'Overwrite?','Yes','No','No');
    if strcmp(answ,'No')
        fprintf(1,'No file created. Have a nice day!\n')
        return
    end
end

csFile = dir(fullfile(expFolder,'*_CondSig.mat'));
if ~isempty(csFile)
    stimSig = load(fullfile(csFile(1).folder,csFile(1).name),'chan*','head*');
else
    try 
        getConditionSignalsBF(fopen(fullfile(smrxFile(1).folder, smrxFile(1).name)));
    catch
        warning('Please create the CondSig.mat file first\n')
        return
    end
end

fsFile = dir(fullfile(expFolder,'*_sampling_frequency.mat'));
if ~isempty(fsFile)
    load(fullfile(fsFile(1).folder,fsFile(1).name),'fs');
else
    warning('Please create the _sampling_frequency.mat file first\n')
    return
end

%%

fields = fieldnames(stimSig);
chanFlag = cellfun(@contains,fields,repmat({'chan'},numel(fields),1));
chanSubs = find(chanFlag);
headers = cellfun(@strrep,fields(chanFlag),...
    repmat({'chan'},numel(chanSubs),1),...
    repmat({'head'},numel(chanSubs),1),'UniformOutput',false);
titles = cell(numel(headers),1);
whiskFlag = false(numel(titles),1);
laserFlag = whiskFlag;
lfpFlag = whiskFlag;
for chead = 1:numel(headers)
    titles{chead} = stimSig.(headers{chead}).title;
    whiskFlag(chead) = checkSignal(titles{chead},'piezo') |...
        checkSignal(titles{chead},'puff');
    laserFlag(chead) = checkSignal(titles{chead},'laser');
    lfpFlag(chead) = checkSignal(titles{chead},'lfp');
end
%% Possible user interaction and correct assignment of the signals
whiskSubs = 1:numel(titles);
laserSubs = whiskSubs;
lfpSubs = whiskSubs;
if any(whiskFlag)
    whiskSubs = find(whiskFlag);
end
while sum(whiskFlag) ~= 1
    wSub = listdlg('ListString',titles(whiskSubs),...
        'PromptString','Select the mechanical TTL',...
        'SelectionMode','single');
    if ~isempty(wSub)
        whiskFlag = false(size(whiskFlag));
        whiskFlag(wSub) = true;
        
    else
        fprintf(1,'Please select one of the displayed signals!\n')
    end
end
whisk = stimSig.(fields{chanSubs(whiskFlag)});
whiskHead = stimSig.(headers{whiskFlag});

if any(laserFlag)
    laserSubs = find(laserFlag);
end
while sum(laserFlag) ~= 1
    lSub = listdlg('ListString',titles(laserSubs),...
        'PromptString','Select the laser TTL',...
        'SelectionMode','single');
    if ~isempty(lSub)
        laserFlag = false(size(laserFlag));
        laserFlag(lSub) = true;
    else
        fprintf(1,'Please select one of the displayed signals!\n')
    end
end
laser = stimSig.(fields{chanSubs(laserFlag)});

if any(lfpFlag)
    lfpSubs = find(lfpFlag);
end
iOk = true;
while sum(lfpFlag) ~= 1 && iOk
    [lSub, iOk] = listdlg('ListString',titles(lfpSubs),...
        'PromptString','Select the lfp signal',...
        'SelectionMode','single','CancelString', 'none');
    if ~isempty(lSub)
        lfpFlag = false(size(lfpFlag));
        lfpFlag(lSub) = true;
    end
end
if iOk
    lfp = stimSig.(fields{chanSubs(lfpFlag)});
else
    lfp = 0;
end

%% Trigger variable construction
intanDomain = round([whiskHead.start, whiskHead.stop]*whiskHead.SamplingFrequency);
if intanDomain(2) <= length(lfp) 
    lfp = lfp(intanDomain(1):intanDomain(2));
elseif diff(intanDomain) <= length(lfp)
    try
        lfp = lfp(intanDomain(1):length(whisk));
    catch
        lfp = lfp(intanDomain(1):length(lfp));
    end
end
Triggers.whisker = whisk; Triggers.laser = laser; Triggers.lfp = lfp;
%% Subscript processing and stimukus finding
wObj = StepWaveform(whisk,fs,'on/off','Mechanical TTL');
lObj = StepWaveform(laser,fs,'on/off','Laser TTL');
wSub = wObj.subTriggers;
lSub = lObj.subTriggers;
glitchInWhisker = diff(wSub,1,2) == 0;
if any(glitchInWhisker)
    fprintf(1,'There were some ''funky'' triggers in whisker... deleting\n')
    wSub(glitchInWhisker,:) = [];
end
glitchInLaser = diff(lSub,1,2) == 0;
if any(glitchInLaser)
    fprintf(1,'There were some ''funky'' triggers in laser... deleting\n')
    lSub(glitchInLaser,:) = [];
end

Conditions(1).name = 'WhiskerAll';
Conditions(1).Triggers = wSub;
Conditions(2).name = 'LaserAll';
Conditions(2).Triggers = lSub;
maxPulses = min(size(lSub,1),size(wSub,1));
dm = distmatrix(lSub(:,1)/fs,wSub(:,1)/fs);
[srtDelay, whr] = sort(dm(:),'ascend');
[~, piezSub] = ind2sub(size(dm),whr(1:maxPulses));
timeDelay = srtDelay(1:maxPulses);
delays = 10.^uniquetol(log10(timeDelay),0.01/log10(max(abs(timeDelay))));
delays(delays > 1) = [];
if std(delays.*1e3) < 1
    delays = mean(delays);
end
Ndel = numel(delays);
fprintf(1,'Delays found:')
lsDel = false(length(timeDelay),Ndel);
Ncond = numel(Conditions);
for cdl = 1:Ndel
    fprintf(1,' %.1f',delays(cdl)*1e3)
    lsDel(:,cdl) = ismembertol(log10(timeDelay),log10(delays(cdl)),...
        abs(0.01/log10(max(delays))));
    Conditions(Ncond + cdl).name = sprintf('Delay %0.3f s',...
        delays(cdl));
    Conditions(Ncond + cdl).Triggers =...
        wSub(sort(piezSub(lsDel(:,cdl))),:);
end
%{
this is a fix to account for a  file with no delays (e.g. period of whisker
stimuli followed by period of laser stimuli)
%}
if Ndel==0,cdl=0;
        fprintf(1,' none!')
else
    fprintf(1,' ms\n')
end 
[~,lghtSub] = min(dm,[],2,'omitnan');
[~,piezSub] = min(dm,[],1,'omitnan');
piezSub = piezSub';
loneLaser = findUnpairedPulse(lghtSub);
lonePiezo = [reshape(diff(sort(piezSub)) == 0,numel(piezSub)-1,1);false];
Conditions(Ncond + cdl + 1).name = 'Laser Control';
Conditions(Ncond + cdl + 1).Triggers = lSub(loneLaser,:);
Conditions(Ncond + cdl + 2).name = 'WhiskerStim Control';
Conditions(Ncond + cdl + 2).Triggers = wSub(lonePiezo,:);
%{
if there are no delays, then WhiskerAll is redundant with WhiskerStim
Control and  LaserAll is redundant with Laser Control, so just take control
condtions
%}
if Ndel==0
    Conditions=Conditions((end-1):end);
end
%%
save(fullfile(expFolder,[expName,'analysis.mat']),'Conditions','Triggers')
end


