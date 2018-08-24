function [tIdx, koIdx] = logicEngine(IDsignal, discreteTraces)
%LOGICENGINE Translates the user defined conditions to a logical train to
%remove an event out of the analysis. To add a time consideration is among
%the next steps. Currently, the function only takes care of the appearence of the
%considered events in the further analyses.
%   The function accepts one input argument: the ID of the signals.
% Emilio Isaias-Camacho GrohLabs 2018

tIdx = selectTrigger(IDsignal);
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
koIdx = false(numel(IDsignal),1);
disp('User prompt: Selection of interesting discrete signals...')
[kickOutIndexes,iok] = listdlg(...
    'PromptString','Select the interesting events:',...
    'ListString',IDsignal,...
    'SelectionMode','multiple',...
    'CancelString','None',...
    'OKString','OK',...
    'Name','Selection of discrete signals',...
    'ListSize',[160,15*numel(IDsignal)]);
if iok
    koIdx(kickOutIndexes) = true;
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

end