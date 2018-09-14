function [tIdx, koIdx, iIdx] = logicEngine(IDsignal, discreteTraces)
%LOGICENGINE Translates the user defined conditions to a logical train to
%select a trigger from the conditioning variables, remove an event out of ,
%or to treat a signal as 'irrelevant' for the analysis.
%   The function accepts one input argument: the ID of the signals, and the
%   discrete time signal array to verify the occurrance of any of the
%   conditioning variables across time.
% Emilio Isaias-Camacho GrohLabs 2018
availableSignals = sum(discreteTraces,2) > 0;
excludeSignal = strcmp(IDsignal,'exclude') | strcmp(IDsignal,'grooming');
TOTAL_EVENTS = numel(IDsignal);
populatedSignals = 1:TOTAL_EVENTS;
populatedSignals(~availableSignals | excludeSignal) = [];
keepSignals = false(TOTAL_EVENTS,1);
triggerSignal = false(TOTAL_EVENTS,1);
ignoreSignal = keepSignals;
tIdx = selectTrigger(IDsignal(populatedSignals));
if ~islogical(tIdx)
    tIdx = -1;
    koIdx = -1;
    return
end
pIdx = checkPossibleCombinations(tIdx,discreteTraces(populatedSignals,:));
chSig = kickOutEvents(IDsignal(populatedSignals(~tIdx & pIdx)),...
    populatedSignals(~tIdx & pIdx));
keepSignals(chSig) = true;
chSig = holyIgnorance(IDsignal(populatedSignals(keepSignals)),chSig);
ignoreSignal(chSig) = true;
triggerSignal(populatedSignals(tIdx)) = true;
if sum(keepSignals & triggerSignal & ignoreSignal)
    error('Something went really wrong. Debugging process required!')
else
    koIdx = keepSignals;
    tIdx = triggerSignal;
    iIdx = ignoreSignal;
end
end

function idx = selectTrigger(IDsignal)
idx = false(numel(IDsignal),1);
[triggerIdx,iok] = listdlg(...
    'PromptString','Select one trigger signal:',...
    'ListString',IDsignal,...
    'SelectionMode','single',...
    'CancelString','None',...
    'OKString','OK',...
    'Name','Selection of discrete signals',...
    'ListSize',[160,15*numel(IDsignal)]);
if iok
    idx(triggerIdx) = true;
else
    idx = -1;
    disp('If you don''t select a trigger, no stack can be built!')
end
end

function chSig = kickOutEvents(IDsignal,inputSignals)
koIdx = true(numel(IDsignal),1);
disp('User prompt: Selection of exclusive signals...')
[kickOutIndexes,iok] = listdlg(...
    'PromptString','Select the excluding signals:',...
    'ListString',IDsignal,...
    'SelectionMode','multiple',...
    'CancelString','None',...
    'OKString','OK',...
    'Name','Selection of discrete signals',...
    'ListSize',[160,15*numel(IDsignal)]);
% What does this mean?
if iok
    koIdx(kickOutIndexes) = false;
else
    
end
chSig = inputSignals(koIdx);
if ~isempty(kickOutIndexes)
    disp('Thank you. You chose the following signals to be excluded:')
    for cid = kickOutIndexes
        fprintf('%s ',IDsignal{cid})
    end
    fprintf('\n');
else
    disp('All available events will be considered for further analysis')
end
end

function chSig = holyIgnorance(IDsignal,signalSub)
iIdx = false(numel(IDsignal,1));
disp('User prompt: Selection of ''irrelevant'' signals:')
[ignoreIndexes, iok] = listdlg(...
    'PromptString','Select the signals to ignore:',...
    'ListString',IDsignal,...
    'SelectionMode','multiple',...
    'CancelString','None',...
    'OKString','OK',...
    'Name','Ignore...',...
    'ListSize',[160,15*numel(IDsignal)]);
if iok
    iIdx(ignoreIndexes) = true;
else
end
chSig = signalSub(iIdx);
if ~isempty(ignoreIndexes)
    disp('Thank you. You chose the following signals to be ignored:')
    for cid = ignoreIndexes
        fprintf('%s ',IDsignal{cid})
    end
    fprintf('\n');
else
    disp('All conditioning variables will be taken into consideration')
    disp('for inclusion or exclusion')
end
end

function pIdx = checkPossibleCombinations(tIdx, discreteTraces)
disp('Exclusive combinations!')
pIdx = ~tIdx;
% evntsPerChannel = sum(discreteTraces(:,discreteTraces(tIdx,:)),2);
evntsPerChannel = sum(discreteTraces,2);
pIdx = pIdx & evntsPerChannel;
end