function kIdx = purgeTrials(discreteStack, timeLapse, tIdx, koIdx, IDs, fs)
%PURGETRIALS Takes the created discrete stack and searches for those trials
%which are not relevant to the experimenter or which need to be excluded to
%proceed with the analysis parting from the stack.
%   
% Emilio Isaias-Camacho @GrohLAb 2018
%% TODO
% 1.- The delayArray should have the data structure such that a 0 delay is
% considered as a synchronous stimuli. The combination with this variable
% and the koIdx variable should be enough to know if a 0 is really a no
% delay or a doesn't matter a delay.
% 2.- This should also be included to the ISI and the triggered average
% computation.

%% FUNCTION
defTL = cell(1,sum(koIdx));
for ctl = 1:sum(koIdx)
    defTL(ctl) = {[num2str(-timeLapse(1)*1e3),', ',num2str(timeLapse(2)*1e3)]};
end
windowTimes = inputdlg(...
    IDs(koIdx),...
    'Delay for the interesting signals',...
    [1,20],...
    defTL);

windowArray = cell2mat(cellfun(@str2num, windowTimes, 'UniformOutput', false));
indexWindow = (windowArray + timeLapse(1)*1e3)*fs*1e-3 + 1;

[Ne, Ns, Na] = size(discreteStack);
if Ne > 2 && sum(~koIdx)
    eOccurrances = squeeze(sum(discreteStack,2));
    kIdx = sum(eOccurrances([false(2,1);~koIdx(~tIdx)],:)) > 0;
    auxIdx = true(1,size(eOccurrances,2));
    auxOcc = eOccurrances([false(2,1);koIdx(~tIdx)],:);
    for ce = 1:sum(koIdx)
        auxIdx = auxIdx & auxOcc(ce,:);
    end
    kIdx = kIdx | ~auxIdx;
end



end