function [eIdxArray, excludeIdx, windowArray] =...
    purgeTrials(discreteStack, timeLapse, tIdx, koIdx, iIdx, IDs, fs)
%PURGETRIALS Takes the created discrete stack and searches for those trials
%which are not relevant to the experimenter or which need to be excluded to
%proceed with the analysis parting from the stack.
%   
% Emilio Isaias-Camacho @GrohLab 2018
%% TODO
% 1.- The delayArray should have the data structure such that a 0 delay is
% considered as a synchronous stimuli. The combination with this variable
% and the koIdx variable should be enough to know if a 0 is really a no
% delay or a doesn't matter a delay.
% 2.- This should also be included to the ISI and the triggered average
% computation.

%% FUNCTION
Na = size(discreteStack,3);
defTL = cell(1,sum(koIdx));
for ctl = 1:sum(koIdx)
    defTL(ctl) = {[num2str(-timeLapse(1)*1e3),', ',num2str(timeLapse(2)*1e3)]};
end

if sum(koIdx & ~iIdx)
    windowTimes = inputdlg(...
        IDs(koIdx & ~iIdx),...
        'Delay for the interesting signals',...
        [1,20],...
        defTL(1:sum(koIdx & ~iIdx)));
    triggerName = IDs{tIdx};
    IDs(tIdx) = [];
    windowArray = cell2mat(cellfun(@str2num, windowTimes, 'UniformOutput', false));
    indexWindow = (windowArray + timeLapse(1)*1e3)*fs*1e-3 + 1;
else
    disp('Either all the signals will be ignored or excluded')
    eIdxArray = koIdx & ~iIdx;
    excludeIdx = ~koIdx;
    windowArray = [1 (sum(timeLapse)*fs + 1)];
    return
end
% Inclusion or exclusion of the signals
eIdxArray = false(sum(koIdx & ~iIdx),Na);
% Event index in the discrete stack
eventIdxDS = find(koIdx(~tIdx) & ~iIdx(~tIdx)) + 2;
excludeIdx =...
    sum(squeeze(sum(discreteStack(...
    [false(2,1);...
    ~koIdx(~tIdx)],...
    :,:),2)),1) > 0;

% For every considered event and every alignment point which meets the
% relaxed criteria we need to verify the windows
deleteEmptyEvents = false(length(eventIdxDS),1);
for ce = 1:length(eventIdxDS)
    eIdxArray(ce,:) =...
        squeeze(sum(discreteStack(...
        eventIdxDS(ce),indexWindow(ce,1):indexWindow(ce,2),:),2)) > 0;
    if ~sum(eIdxArray(ce,:))
        fprintf('%s has no action between %.3f and %.3f ms relative to %s\n',...
            IDs{eventIdxDS(ce)-2},windowArray(ce,1),windowArray(ce,2),...
            triggerName)
        disp('Deleting...')
        deleteEmptyEvents(ce) = true;
    end
end
eIdxArray(deleteEmptyEvents,:) = [];

end