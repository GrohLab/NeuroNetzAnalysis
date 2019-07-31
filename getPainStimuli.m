function [Conditions, Triggers] = getPainStimuli(expFolder)
% GETPAINSTIMULI searches specifically for the stimulation conditions for
% the CFA/Saline mechanical pressure and laser on the somatosensory cortex
% of the hindlimb. The existing conditions are mechanical pressure plus 1
% Hz laser, mechanical pressure plus 10 Hz laser, mechanical control and
% laser control. 
checkSignal = @(x,y) contains(x,y,'IgnoreCase',true);
Conditions = struct('name',{},'Triggers',{});
Triggers = struct();
csFile = dir(fullfile(expFolder,'*_CondSig.mat'));
if ~isempty(csFile)
    stimSig = load(fullfile(csFile.folder,csFile.name),'chan*','head*');
else
    warning('Please create the CondSig.mat file first\n')
    return
end

fsFile = dir(fullfile(expFolder,'*_sampling_frequency.mat'));
if ~isempty(fsFile)
    load(fullfile(fsFile.folder,fsFile.name),'fs');
else
    warning('Please create the _sampling_frequency.mat file first\n')
    return
end

binFile = dir(fullfile(expFolder,'*.smrx'));
if isempty(binFile)
    warning('There is no experiment in this folder...\n')
    return
end

fields = fieldnames(stimSig);
chanFlag = cellfun(@contains,fields,repmat({'chan'},numel(fields),1));
chanSubs = find(chanFlag);
headers = cellfun(@strrep,fields(chanFlag),...
    repmat({'chan'},numel(chanSubs),1),...
    repmat({'head'},numel(chanSubs),1),'UniformOutput',false);
titles = cell(numel(headers),1);
mechFlag = false(numel(titles),1);
laserFlag = mechFlag;
for chead = 1:numel(headers)
    titles{chead} = stimSig.(headers{chead}).title;
    mechFlag(chead) = checkSignal(titles{chead},'mech');
    laserFlag(chead) = checkSignal(titles{chead},'laser');
end
%% Possible user interaction
mechSubs = find(mechFlag);
while sum(mechFlag) > 1
    mSub = listdlg('ListString',titles(mechFlag),...
        'PromptString','Select the mechanical TTL',...
        'SelectionMode','single');
    if ~isempty(mSub)
        mechFlag = false(size(mechFlag));
        mechFlag(mechSubs(mSub)) = true;
    else
        fprintf(1,'Please select one of the displayed signals!\n')
    end
end
mech = stimSig.(fields{chanSubs(mechFlag)});

laserSubs = find(laserFlag);
while sum(laserFlag) > 1
    lSub = listdlg('ListString',titles(laserFlag),...
        'PromptString','Select the laser TTL',...
        'SelectionMode','single');
    if ~isempty(lSub)
        laserFlag = false(size(laserFlag));
        laserFlag(laserSubs(lSub)) = true;
    else
        fprintf(1,'Please select one of the displayed signals!\n')
    end
end
laser = stimSig.(fields{chanSubs(laserFlag)});
%% Subscript processing
mObj = StepWaveform(mech,fs,'on/off','Mechanical TTL');
lObj = StepWaveform(laser,fs,'on/off','Laser TTL');
mSub = mObj.subTriggers;
lSub = lObj.subTriggers;

lsIpi = diff(lSub(:,1)/fs);
lsFst = StepWaveform.firstOfTrain(lSub(:,1)/fs, 5 - 1e-3);
lsIpi = lsIpi(~lsFst(1:end-1));
pulsFreq = 1./diff(lSub(:,1)./fs);
freqCond = round(uniquetol(pulsFreq,0.1/max(pulsFreq)));
freqCond = unique(freqCond(freqCond > 0));
fprintf(1,'Frequency stimulation:')

if isempty(freqCond) || numel(freqCond) == 1
    freqCond = 0;
    fprintf(' None')
elseif numel(freqCond) > 1
    Nfre = numel(freqCond);
    lsuSub = lSub(lsFst,1);
    lsdSub = lSub(:,2);
    lsLst = flip(StepWaveform.firstOfTrain(flip(lsdSub)/fs, 5-1e-3));
    lsdSub = lsdSub(lsLst);
    lsCon = [lsuSub,lsdSub];
    lsIdx = false(size(lsCon,1),Nfre);
    pulseFreqTrain = pulsFreq(lsFst);
    NT = 0;
    for cdl = 1:Nfre
        fprintf(1,' %.2f',freqCond(cdl))
        lsIdx(:,cdl) = ismembertol(pulseFreqTrain,freqCond(cdl),...
            0.01*max(pulseFreqTrain));
        Conditions(cdl).name = sprintf('%.2f Hz',freqCond(cdl));
        isd = distmatrix(lsCon(lsIdx(:,cdl),1),mSub(:,1));
        [pSub, ~] = find(isd < fs*0.0005);
        plSub = lsCon(lsIdx(:,cdl),:);
        Conditions(cdl).Triggers = plSub(pSub,:);
        NT = NT + size(plSub(pSub,:),1);
    end
end
fprintf(1,' Hz\n')
growingSubs = zeros(NT,2);
k = 0;
for cpc = 1:numel(Conditions)
    csz = size(Conditions(cpc).Triggers,1);
    growingSubs(k+1:k+csz,:) = Conditions(cpc).Triggers;
    k = csz;
end
NpairCond = numel(Conditions);
Conditions(NpairCond+1).name = 'Laser Control';
Conditions(NpairCond+1).Triggers = setdiff( lsCon, growingSubs,'rows');

isd = distmatrix(mSub(:,1), sort(growingSubs(:,1)));
[pSub, ~] = find(isd < fs*0.0005);
pmFlag = true(size(mSub,1),1);
pmFlag(pSub) = false;
Conditions(NpairCond+2).name = 'Mech Control';
Conditions(NpairCond+2).Triggers = mSub(pmFlag,:);

Triggers = struct('whisker',mech,'laser',laser);
[~,expName,~]=fileparts(binFile(1).name);
save(fullfile(expFolder,[expName,'analysis.mat']),'Conditions','Triggers')
end


