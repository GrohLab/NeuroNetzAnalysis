%% Experiment pooling
% Experiment folder search
experimentDir = uigetdir('Z:\Jesus\LTP_Jesus_Emilio',...
    'Select an experiment directory');
% Folders existing in the selected directory
expFolders = dir(experimentDir);
expFolders(1:2) = [];
% Selection of the considered experiments
[chExp, iOk] = listdlg('ListString', arrayfun(@(x) x.name, expFolders,...
    'UniformOutput', 0), 'InitialValue', 1:size(expFolders,1),...
    'PromptString', 'Select the experiments to consider:',...
    'SelectionMode', 'multiple', 'ListSize', [350, (numel(expFolders))*16]);
Nexp = numel(chExp);

%% User controlling variables
% Time lapse, bin size, and spontaneous and response windows
promptStrings = {'Viewing window (time lapse) [s]:','Response window [s]',...
    'Bin size [s]:'};
defInputs = {'-0.1, 0.1', '0.002, 0.05', '0.001'};
answ = inputdlg(promptStrings,'Inputs', [1, 30],defInputs);
if isempty(answ)
    fprintf(1,'Cancelling...\n')
    return
else
    timeLapse = str2num(answ{1}); %#ok<*ST2NM>
    if numel(timeLapse) ~= 2
        timeLapse = str2num(inputdlg('Please provide the time window [s]:',...
            'Time window',[1, 30], '-0.1, 0.1'));
        if isnan(timeLapse) || isempty(timeLapse)
            fprintf(1,'Cancelling...')
            return
        end
    end
    responseWindow = str2num(answ{2});
    binSz = str2double(answ(3));
end
fprintf(1,'Time window: %.2f - %.2f ms\n',timeLapse(1)*1e3, timeLapse(2)*1e3)
fprintf(1,'Response window: %.2f - %.2f ms\n',responseWindow(1)*1e3, responseWindow(2)*1e3)
fprintf(1,'Bin size: %.3f ms\n', binSz*1e3)
sponAns = questdlg('Mirror the spontaneous window?','Spontaneous window',...
    'Yes','No','Yes');
spontaneousWindow = -flip(responseWindow);
if strcmpi(sponAns,'No')
    spontPrompt = "Time before the trigger in [s] (e.g. -0.8, -0.6 s)";
    sponDef = string(sprintf('%.3f, %.3f',spontaneousWindow(1),...
        spontaneousWindow(2)));
    sponStr = inputdlg(spontPrompt, 'Inputs',[1,30],sponDef);
    if ~isempty(sponStr)
        spontAux = str2num(sponStr{1});
        if length(spontAux) ~= 2 || spontAux(1) > spontAux(2) || ...
                spontAux(1) < timeLapse(1)
            fprintf(1, 'The given input was not valid.\n')
            fprintf(1, 'Keeping the mirror version!\n')
        else
            spontaneousWindow = spontAux;
        end
    end
end
fprintf(1,'Spontaneous window: %.2f to %.2f ms before the trigger\n',...
    spontaneousWindow(1)*1e3, spontaneousWindow(2)*1e3)
statFigFileNameEndings = {'.pdf','.emf'};
figFormat = {'-dpdf','-dmeta'};

cellLogicalIndexing = @(x,idx) x(idx);

expCo = 1;
%% Experiment loop 
% For all the chosen experiments, perform the statistical test and save its
% results, and extract the spike times for all clusters blending in the
% trial information for a PDF calculation
relativeSpikeTimes = cell(Nexp);
for cexp = chExp
    if ~loadTriggerData(expFolders(cexp).name)
        fprintf(1,'Skipping experiment %s\n', expFolders(cexp).name)
        continue
    end
    % Figures directory for the considered experiment
    figureDir = fullfile(dataDir,'Figures\');
    if ~mkdir(figureDir)
        fprintf(1,'There was an issue with the figure folder...\n');
    end
    Ns = min(structfun(@numel,Triggers));
    % Total duration of the recording
    Nt = Ns/fs;
    % Useless clusters (labeled as noise or they have very low firing rate)
    badsIdx = cellfun(@(x) x==3,sortedData(:,3));
    bads = find(badsIdx);
    totSpkCount = cellfun(@numel,sortedData(:,2));
    clusterSpikeRate = totSpkCount/Nt;
    silentUnits = clusterSpikeRate < 0.1;
    bads = union(bads,find(silentUnits));
    goods = setdiff(1:size(sortedData,1),bads);
    badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));
    % Logical spike trace for the first good cluster
    spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
    % Subscript column vectors for the rest good clusters
    spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods(2:end),2),...
        'UniformOutput',false);
    % Number of good clusters
    Ncl = numel(goods);
    gclID = sortedData(goods,1);
    continuousSignals = {Triggers.whisker; Triggers.laser};
    % User defined conditions variables if not already chosen
    if ~exist('alignCond','var') && ~exist('Conditions','var')
        condNames = arrayfun(@(x) x.name, Conditions, 'UniformOutput', false);
        condGuess = contains(condNames, 'whiskerall', 'IgnoreCase', true);
        % Choose the conditions to create the stack upon
        [alignCond, iOk] = listdlg('ListString',condNames,'SelectionMode','single',...
            'PromptString',...
            'Choose the condition which has all whisker triggers: (one condition)',...
            'InitialValue', find(condGuess), 'ListSize', [350, numel(condNames)*16]);
        if ~iOk
            fprintf(1,'Cancelling...\n')
            return
        end
        % Select the onset or the offset of a trigger
        fprintf(1,'Condition ''%s''\n', Conditions(alignCond).name)
        onOffStr = questdlg('Trigger on the onset or on the offset?','Onset/Offset',...
            'on','off','Cancel','on');
        if strcmpi(onOffStr,'Cancel')
            fprintf(1,'Cancelling...\n')
            return
        end

        % Choose the conditions to look at
        auxSubs = setdiff(1:numel(condNames), alignCond);
        ccondNames = condNames(auxSubs);
        [cchCond, iOk] = listdlg('ListString',ccondNames,'SelectionMode','multiple',...
            'PromptString',...
            'Choose the condition(s) to look at (including whiskers):',...
            'ListSize', [350, numel(condNames)*16]);
        if ~iOk
            fprintf(1,'Cancelling...\n')
            return
        end
        
        % Select the onset or the offset of a trigger
        fprintf(1,'Condition(s):\n')
        fprintf('- ''%s''\n', Conditions(auxSubs(cchCond)).name)
        % Subscripts and names for the considered conditions
        consideredConditions = auxSubs(cchCond);
        consCondNames = arrayfun(@(x) x.name, Conditions(consideredConditions),...
            'UniformOutput', 0);
        Nccond = length(consideredConditions);
        
        sponActStackIdx = tx >= spontaneousWindow(1) & tx <= spontaneousWindow(2);
        respActStackIdx = tx >= responseWindow(1) & tx <= responseWindow(2);
        % The spontaneous activity of all the clusters, which are allocated from
        % the second until one before the last row, during the defined spontaneous
        % time window, and the whisker control condition.
        
        timeFlags = [sponActStackIdx;respActStackIdx];
        % Time window
        delta_t = diff(responseWindow);
        % Spontaneous vs evoked comparison
        indCondSubs = cumsum(Nccond:-1:1);
        isWithinResponsiveWindow =...
            @(x) x > responseWindow(1) & x < responseWindow(2);
    end
    % Constructing the stack out of the user's choice
    % discStack - dicrete stack has a logical nature
    % cst - continuous stack has a numerical nature
    % Both of these stacks have the same number of time samples and trigger
    % points. They differ only in the number of considered events.
    [discStack, cst] = getStacks(spkLog,Conditions(alignCond).Triggers,onOffStr,...
        timeLapse,fs,fs,spkSubs,continuousSignals);
    % Number of clusters + the piezo as the first event + the laser as the last
    % event, number of time samples in between the time window, and number of
    % total triggers.
    [Ne, Nt, NTa] = size(discStack);
    % Computing the time axis for the stack
    tx = (0:Nt - 1)/fs - timeLapse(1);
    
    % Computing which alignment points belong to which condition.
    delayFlags = false(NTa,Nccond);
    counter2 = 1;
    for ccond = consideredConditions
        delayFlags(:,counter2) = ismember(Conditions(alignCond).Triggers(:,1),...
            Conditions(ccond).Triggers(:,1));
        counter2 = counter2 + 1;
    end
    Na = sum(delayFlags,1);

    % Statistical tests
    [Results, Counts] = statTests(discStack, delayFlags, timeFlags);
    
    % Plotting statistical tests
    figs = scatterSignificance(Results, Counts, consCondNames, delta_t, gclID);
    arrayfun(@configureFigureToPDF, figs);
    % Firing rate for all clusters, for all trials
    % meanfr = cellfun(@(x) mean(x,2)/delta_t,Counts,'UniformOutput',false);
    for cfig = 1:numel(figs)
        statFigName = strsplit(figs(cfig).Children(1).Title.String, ': ');
        statFigName = statFigName{2};
        statFigName = [statFigName,...
            sprintf(' RW:%.1f-%.1f ms', responseWindow(1)*1e3, responseWindow(2)*1e3)];
        print(figs(cfig), fullfile(figureDir, [statFigName, '.pdf']),...
            '-dpdf','-fillpage');
        print(figs(cfig), fullfile(figureDir, [statFigName, '.emf']),...
            '-dmeta');
        close(figs(cfig))
    end
    
    % Filtering for the whisker responding clusters
    H = cell2mat(cellfun(@(x) x.Pvalues,...
        arrayfun(@(x) x.Activity, Results(indCondSubs), 'UniformOutput', 0),...
        'UniformOutput', 0)) < 0.05;
    Htc = sum(H,2);
    % Those clusters responding more than 80% of all whisker stimulating
    % conditions
    wruIdx = Htc/Nccond > 0.80;
    Nwru = nnz(wruIdx);
    gclID = sortedData(goods,1);
    fprintf('%d whisker responding clusters:\n', Nwru);
    fprintf('- %s\n',gclID{wruIdx})
    if ~Nwru
        cans = questdlg('No evoked response! Skip?','No response!',...
            'Yes','No','Yes');
        if strcmp(cans,'Yes')
            continue
        end
    end
    
    % Getting spike times for every responding cluster and computing its
    % PDF
    
    for ccond = 1:size(delayFlags,2)
        relativeSpikeTimes(expCo) = getRasterFromStack(discStack,~delayFlags(:,ccond),...
            wruIdx, timeLapse, fs, true, false);
        relativeSpikeTimes{expCo}(:,~delayFlags(:,ccond)) = [];
        respIdx = cellfun(isWithinResponsiveWindow, relativeSpikeTimes{expCo},...
            'UniformOutput', false);
        spikeTimesINRespWin = cellfun(cellLogicalIndexing,...
            relativeSpikeTimes{expCo}, respIdx, 'UniformOutput',false);
        
        
    end
    expCo = expCo + 1;
    
end