function koIdx = logicEngine(discreteStackSize,IDsignal)
%LOGICENGINE Translates the user defined conditions to a logical train to
%remove an event out of the analysis. To add a time consideration is among
%the next steps. Currently, the function only takes care of the appearence of the
%considered events in the further analyses.
%   The function accepts two input arguments. The discrete stack size and
%   the ID of the signals. The size is only to add a verification step that
%   would save some trouble.
if discreteStackSize(1) - 2 == numel(IDsignal)
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