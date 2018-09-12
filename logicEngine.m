function [tIdx, koIdx] = logicEngine(IDsignal, discreteTraces, fs)
%LOGICENGINE Translates the user defined conditions to a logical train to
%remove an event out of the analysis. To add a time consideration is among
%the next steps. Currently, the function only takes care of the appearence of the
%considered events in the further analyses.
%   The function accepts one input argument: the ID of the signals.
% Emilio Isaias-Camacho GrohLabs 2018

TOTAL_EVENTS = numel(IDsignal);
availableTrig = true(TOTAL_EVENTS,1);
tp = 0;
populatedSignals = 1:TOTAL_EVENTS;
while sum(availableTrig) <= TOTAL_EVENTS && sum(availableTrig) > 0 &&...
        tp == 0
    tIdx = selectTrigger(IDsignal(availableTrig));
    if islogical(tIdx)
        stWf = StepWaveform(discreteTraces(tIdx,:), fs);
        tp = numel(stWf.Triggers);
        if ~tp
            populatedSignals(tIdx) = [];
            discreteTraces(tIdx,:) = [];
            IDsignal(tIdx) = [];
            AVAIL_EVENTS = numel(IDsignal);
            availableTrig = true(AVAIL_EVENTS,1);
        end
    else
        tIdx = -1;
        koIdx = -1;
        return
    end
end
pIdx = checkPossibleCombinations(tIdx,discreteTraces);
koIdx = kickOutEvents(IDsignal(~tIdx & pIdx));
koSubs = find(koIdx);
koSubs(koSubs >= find(tIdx)) = koSubs(koSubs >= find(tIdx)) + 1;
koIdx = false(numel(IDsignal),1);
koIdx(koSubs) = true;

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

function koIdx = kickOutEvents(IDsignal)
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
if iok
    koIdx(kickOutIndexes) = false;
else
    
end
if ~isempty(kickOutIndexes)
    disp('Thank you. You chose the following signals as interesting:')
    for cid = kickOutIndexes
        fprintf('%s ',IDsignal{cid})
    end
    fprintf('\n');
else
    disp('No other event will be considered for further analysis')
end
end

function pIdx = checkPossibleCombinations(tIdx, discreteTraces)
disp('Exclusive combinations!')
pIdx = ~tIdx;
evntsPerChannel = sum(discreteTraces(:,discreteTraces(tIdx,:)),2);
pIdx = pIdx & evntsPerChannel;
end