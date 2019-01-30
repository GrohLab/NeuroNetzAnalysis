function [configStruct] = createConfigStruct(EphysDir)
%CREATECOFIGSTRUCT returns a structure with the information required to
%proceed with a population analysis.
%   Prompting the user for specific inputs. For example, the cell type can
%   be red from the RecDB variable. However, the conditional variables or
%   signals need to be assumed and hard-coded in this function. 
dbTabFile = fullfile(EphysDir,'ephys_database.mat');
deleteLastCell = @(x) x(1:end-1);
if exist(dbTabFile,'file')
    anMatDir = fullfile(EphysDir,'EphysData','AnalysisMatFiles');
    load(dbTabFile,'RecDB')
    Nx = size(RecDB,1);
    cellTypes = unique(RecDB.PhysioNucleus);
    %% Cell type selection
    fFlag = false(Nx,1);
    [answ, lIdx] = listdlg('ListString',cellTypes,'PromptString',...
        'Select a cell type:','SelectionMode','multiple',...
        'ListSize',[160,15*numel(cellTypes)],'CancelString','None',...
        'OKString','OK','Name','Cell type selection');
    if ~lIdx
        [sfans] = questdlg('Would you like to perform single file analysis?',...
            'Alternative analysis','Yes','No','Yes');
        if ~strcmpi(sfans,'yes')
            configStruct = [];
            return
        end
        [fName, fDir] = uigetfile('*analysis.mat',...
            'Select an analysis file',anMatDir);
        % Give the animal name and the start time of the selected
        % experiment
        expNames = deleteLastCell(strsplit(fName,'analysis.mat'));
        fFlag(strcmp(RecDB.Properties.RowNames,expNames)) = true;
        expFileNames = fullfile(fDir,fName);
    end
    % Get the file names from the RecDB or from the single file selector.
    if ~isempty(answ)
        auxCT = cell(numel(answ),1);
        for cct = answ
            fFlag(RecDB.PhysioNucleus == cellTypes(cct)) = true;
            auxCT(cct) = {cellTypes(cct)};
        end
        configStruct.CellTypes = auxCT;
        expNames = RecDB.Properties.RowNames(fFlag);
        % Experiment names to process:
        expFileNames = cellfun(@fullfile,...
            repmat({anMatDir},sum(fFlag),1),...
            cellfun(@strcat,expNames,...
            repmat({'analysis.mat'},sum(fFlag),1), 'UniformOutput',false),...
            'UniformOutput',false);
    end
    Nf = sum(fFlag);
    %% Viewing window definition
    errorTL = true;
    while errorTL
        promtMsg =...
            {'Time Before the trigger, time After the trigger [ms]:'};
        titleMsg = 'Viewing window';
        defInput = {'250, 500'};
        opts = struct('Resize','on','WindowStyle','modal');
        vwAns = inputdlg(promtMsg,titleMsg,[1, 29],defInput,opts);
        timeLapse = str2double(strsplit(vwAns{:},{',',' '}));
        if numel(timeLapse) == 2 
            errorTL = false;
        elseif isempty(vwAns)
            timeLapse = [];
        else
            continue
        end
        break
    end
    % Conversion to seconds
    configStruct.ViewWindow = timeLapse*1e-3;
    %% Bin size definition
    configStruct.BinSize = [];
    promtMsg = 'Bin size [ms]:';
    titleMsg = 'Bin size';
    defInput = {'10'};
    bzAns = inputdlg(promtMsg,titleMsg,[1, 15],defInput,opts);
    bnSz = str2double(bzAns);
    % Conversion to seconds
    if ~isempty(bnSz) && ~isnan(bnSz)
        configStruct.BinSize = bnSz*1e-3;
    end
    %% Trigger selection
    % Hard coded conditioning variables. Future modification for
    % data-dependent names.----------------------------------------------!!
    % 'exclude' and 'grooming' are not presented to the user by default.
    condVars = {'whisker';'light';'puff';'touch';'pole'};
    [triggerIdx,iok] = listdlg(...
        'PromptString','Select one trigger signal:',...
        'ListString',condVars,...
        'SelectionMode','single',...
        'CancelString','None',...
        'OKString','OK',...
        'Name','Selection of discrete signals',...
        'ListSize',[160,15*numel(condVars)],'InitialValue',1);
    if iok
        edAns = questdlg('On which edge should the trigger be considered?',...
            'Trigger edge','On-set','Off-set','On-set');
        if strcmp(edAns,'On-set') || isempty(edAns)
            edBool = true;
        else
            edBool = false;
        end
        configStruct.Trigger = struct('Name',condVars(triggerIdx),...
            'Edge',edBool);
    else
        % No trigger selected
        warning('No trigger selected! Exiting.')
        configStruct = [];
        return
    end
    %% Exclude selection
    tIdx = cellfun(@isequal,strfind(condVars,condVars(triggerIdx)),...
        repmat({1},numel(condVars),1));
    exclVars = condVars(~tIdx);
    [exclSub, iok] = listdlg('ListString',exclVars,...
        'PromptString','Select those signals to exclude:',...
        'ListString',exclVars,...
        'SelectionMode','multiple',...
        'CancelString','None',...
        'OKString','OK',...
        'Name','Selection of exclude signals',...
        'ListSize',[160,15*numel(exclVars)],'InitialValue',1);
    if iok
        configStruct.Exclude = exclVars(exclSub);
    else
        fprintf('No signal will be excluded from the experiment')
        configStruct.Exclude = [];
    end
    %% Ignore selection
    eIdx = false(numel(exclVars),1);
    eIdx(exclSub) = true;
    ignVars = exclVars(~eIdx);
    configStruct.Ignore = [];
    if ~isempty(ignVars)
        [ignSub, iok] = listdlg('ListString',ignVars,...
            'PromptString','Select those signals to ignore:',...
            'ListString',ignVars,...
            'SelectionMode','multiple',...
            'CancelString','None',...
            'OKString','OK',...
            'Name','Selection of exclude signals',...
            'ListSize',[160,15*numel(exclVars)],'InitialValue',1);
        if iok
            configStruct.Ignore = ignVars(ignSub);
        end
    end
    %% Conditioning Window selection
    % Include those signals which were not excluded, not ignored, and not
    % being the trigger. Obviously
    ignVars(ignSub)
    exclVars(exclSub)
    condVars(triggerIdx)
else
    
end

end