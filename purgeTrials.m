function [eIdxArray, excludeIdx, windowArray] =...
    purgeTrials(...
    discreteStack,... boolean stack containing all the binary information relative to a trigger
    timeLapse,...     size for the viewing window
    tIdx,...          trigger selection index
    koIdx,...         conditioning variables index
    iIdx,...          non-affecting variables (don't cares)
    IDs,...           names for all the variables in the experiment
    fs...            sampling frequency
    )

%PURGETRIALS Takes the created discrete stack and searches for those trials
%which are not relevant to the experimenter or which need to be excluded to
%proceed with the analysis parting from the stack. It accepts the boolean
%stack, _discreteStack_; the viewing window array with the time before the
%trigger as the first element and the time after the trigger as the second
%element, _timeLapse_; the three boolean indices for trigger selection,
%_tIdx_, for variable exclusion, _koIdx_, and for non-affecting variables,
%_iIdx_, which the further analysis will nor exclude nor actively include;
%the sampling frequency, _fs_; and the _userFlag_, which can be omitted
%from the argument calling. If this last variable is true, the user will be
%able to select the windows for each conditioning variable. Otherwise, the
%function will decide to select half of the viewing window as the
%conditioning window i.e. half time before, half time after.
%
%An additional functionality should be implemented. The ability to load the
%conditioning windows for the non-excluded variables from the conditioning
%window. Maybe we can use a very simple bypass for all the functions.
%
% Emilio Isaias-Camacho @GrohLab 2018

%% Window variable initialization
% The default durations before and after the trigger point are set to the
% complete viewing window.
Na = size(discreteStack,3);
defTL = cell(1,sum(koIdx));
excludeIdx =...
    sum(squeeze(sum(discreteStack(...
    [false(2,1);...
    ~koIdx(~tIdx)],...
    :,:),2)),1) > 0;
for ctl = 1:sum(koIdx)
    defTL(ctl) = {[num2str(-timeLapse(1)*1e3),', ',num2str(timeLapse(2)*1e3)]};
end
%% Window determination
% If the function is called for the user to choose the relative durations
% surrounding the trigger point, a prompting window asks the user for
% specific window durations. In the other case, when the function is called
% to create the 'Control conditions', the functions assumes a conditioning
% window of 50% the size of the 'Viewing window'.
if sum(koIdx & ~iIdx)
    windowTimes = inputdlg(...
        IDs(koIdx & ~iIdx),...
        'Delay for the interesting signals',...
        [1,20],...
        defTL(1:sum(koIdx & ~iIdx)));
    triggerName = IDs{tIdx};
    IDs(tIdx) = [];
    windowArray =...
        cell2mat(...
        cellfun(@str2num, windowTimes, 'UniformOutput', false));
    %         indexWindow = (windowArray + timeLapse(1)*1e3)*fs*1e-3 + 1;
else
    disp('Either all the signals will be ignored or excluded')
    eIdxArray = false(1,Na);
    excludeIdx =...
        sum(squeeze(sum(discreteStack(...
        [false(2,1);...
        ~koIdx(~tIdx)],...
        :,:),2)),1) > 0;
    windowArray = [1 (sum(timeLapse)*fs + 1)];
    return
end
%% Index variables initialization
% Indices for the time axis of the stack
indexWindow = (windowArray*1e-3 + timeLapse(1))*fs + 1;
% Inclusion or exclusion of the signals
eIdxArray = false(sum(koIdx & ~iIdx),Na);
% Event index in the discrete stack
eventIdxDS = find(koIdx(~tIdx) & ~iIdx(~tIdx)) + 2;
deleteEmptyEvents = false(length(eventIdxDS),1);

%% Binary labelling of the triggered trials
% For every considered event and every alignment point which meets the
% relaxed criteria we need to verify the windows
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