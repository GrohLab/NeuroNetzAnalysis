clearvars
%% Experiment pooling
% Experiment folder search
experimentDir = uigetdir('Z:\Jesus\LTP_Jesus_Emilio',...
    'Select an experiment directory');
if experimentDir == 0
    return
end
% Folders existing in the selected directory
expFolders = dir(experimentDir);
expFolders(1:2) = [];
folderNames = arrayfun(@(x) x.name, expFolders, 'UniformOutput', 0);
folderFlag = arrayfun(@(x) x.isdir, expFolders);
% Selection of the considered experiments
[fSel, iOk] = listdlg('ListString', folderNames(folderFlag),...
    'InitialValue', 1:nnz(folderFlag),...
    'PromptString', 'Select the experiments to consider:',...
    'SelectionMode', 'multiple', 'ListSize', [350, (nnz(folderFlag))*16]);
if ~iOk
    fprintf(1, 'Cancelling...\n');
    return
end
folderSub = find(folderFlag);
chExp = folderSub(fSel);
Nexp = numel(chExp);
% Auxiliary variable for the phy multiple updates
phyNames = ["ch","channel";...
    "sh", "shank";...
    "fr", "firing_rate"];

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
printOpts = {{'-dpdf','-fillpage'},'-dmeta'};

cellLogicalIndexing = @(x,idx) x(idx);

expCo = 1;
%% Experiment loop 
% For all the chosen experiments, perform the statistical test and save its
% results, and extract the spike times for all clusters blending in the
% trial information for a PDF calculation
relativeSpikeTimes = cell(Nexp);
for cexp = reshape(chExp, 1, [])
    dataDir = fullfile(expFolders(cexp).folder, expFolders(cexp).name);
    foldContents = dir(dataDir); foldContents(1:2) = [];
    contNames = arrayfun(@(x) x.name, foldContents, 'UniformOutput', 0);
    [~,~,fext] = cellfun(@fileparts, contNames, 'UniformOutput', 0);
    ksPhyFlags = cell2mat(cellfun(@(x) strcmpi(x,{'.npy','.tsv','.py'}),...
        fext, 'UniformOutput', 0));
    dirNames = contNames([foldContents.isdir]);
    if ~any(ksPhyFlags(:))
        if ~isempty(dirNames)
            [subFoldSel, iOk] = listdlg(...
                'ListString', dirNames,...
                'SelectionMode', 'single',...
                'Name', 'Select the experiment folder:',...
                'CancelString', 'Skip');
            if ~iOk
                fprintf(1, 'Skipping experiment %s', expFolders(cexp).name);
                continue
            end
            dataDir = fullfile(dataDir, dirNames{subFoldSel});
        else
            fprintf(1,'This experiment hasn''t been spike-sorted!\n')
            fprintf(1,'Skipping experiment %s\n', expFolders(cexp).name)
            continue
        end
    end
    if ~loadTriggerData(dataDir)
        fprintf(1,'Skipping experiment %s\n', expFolders(cexp).name)
        continue
    end
    %{
    % Figures directory for the considered experiment
    figureDir = fullfile(dataDir,'Figures\');
    if ~mkdir(figureDir)
        fprintf(1,'There was an issue with the figure folder...\n');
        fprintf(1,'Saving the figures in the data folder!\n');
        fprintf(1,'%s\n',dataDir);
        figureDir = dataDir;
    end
    %}
    %% Constructing the helper 'global' variables
    % Number of total samples
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
    badsIdx = badsIdx | silentUnits;
    if ~any(ismember(clInfo.Properties.VariableNames,'ActiveUnit'))
        try
            clInfo = addvars(clInfo,~badsIdx,'After','id',...
                'NewVariableNames','ActiveUnit');
        catch
            clInfo = addvars(clInfo,false(size(clInfo,1),1),'After','id',...
                'NewVariableNames','ActiveUnit');
            clInfo{sortedData(~badsIdx,1),'ActiveUnit'} = true;
            fprintf(1,'Not all clusters are curated!\n')
            fprintf(1,'%s\n',dataDir)
        end
        try
            writeClusterInfo(clInfo,fullfile(dataDir,'cluster_info.tsv'),true);
        catch
            fprintf(1,'Unable to write cluster info for %s\n',dataDir)
        end
    end
    gclID = sortedData(goods,1);
    badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));
    % Logical spike trace for the first good cluster
    spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
    % Subscript column vectors for the rest good clusters
    spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods(2:end),2),...
        'UniformOutput',false);
    % Number of good clusters
    Ncl = numel(goods);
    % Redefining the stimulus signals from the low amplitude to logical values
    whStim = {'piezo','whisker','mech','audio'};
    cxStim = {'laser','light'};
    lfpRec = {'lfp','s1','cortex','s1lfp'};
    trigNames = fieldnames(Triggers);
    numTrigNames = numel(trigNames);
    continuousSignals = struct2cell(Triggers);
    % User defined conditions variables if not already chosen
    if ~exist('chCond','var')
        Nt = round(sum(ceil(abs(timeLapse)*fs))+1);
        % Computing the time axis for the stack
        tx = (0:Nt - 1)/fs + timeLapse(1);
        %% Condition triggered stacks
        condNames = arrayfun(@(x) x.name,Conditions,'UniformOutput',false);
        condGuess = contains(condNames, 'whiskerall', 'IgnoreCase', true);
        % Choose the conditions to create the stack upon
        [chCond, iOk] = listdlg('ListString',condNames,'SelectionMode','single',...
            'PromptString',...
            'Choose the condition which has all whisker triggers: (one condition)',...
            'InitialValue', find(condGuess), 'ListSize', [350, numel(condNames)*16]);
        if ~iOk
            fprintf(1,'Cancelling...\n')
            return
        end
        
        % Select the onset or the offset of a trigger
        fprintf(1,'Condition ''%s''\n', Conditions(chCond).name)
        onOffStr = questdlg('Trigger on the onset or on the offset?','Onset/Offset',...
            'on','off','Cancel','on');
        if strcmpi(onOffStr,'Cancel')
            fprintf(1,'Cancelling...\n')
            return
        end
        
        %% Considered conditions selection
        % Choose the conditions to look at
        auxSubs = setdiff(1:numel(condNames), chCond);
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
        
        ansFilt = questdlg('Would you like to filter for significance?','Filter',...
            'Yes','No','Yes');
        filtStr = 'unfiltered';
        if strcmp(ansFilt,'Yes')
            filtStr = 'filtered';
        end
        % Subscript to indicate the conditions with all whisker stimulations,
        % whisker control, laser control, and the combination whisker and laser.
        consideredConditions = auxSubs(cchCond);
        Nccond = length(consideredConditions);
        
        % Select the onset or the offset of a trigger
        fprintf(1,'Condition(s):\n')
        fprintf('- ''%s''\n', Conditions(auxSubs(cchCond)).name)
        % Subscripts and names for the considered conditions
        consCondNames = condNames(consideredConditions);
        
        % Time windows for comparison between conditions and activity
        sponActStackIdx = tx >= spontaneousWindow(1) & tx <= spontaneousWindow(2);
        respActStackIdx = tx >= responseWindow(1) & tx <= responseWindow(2);
        % The spontaneous activity of all the clusters, which are allocated from
        % the second until one before the last row, during the defined spontaneous
        % time window, and the whisker control condition.
        
        timeFlags = [sponActStackIdx; respActStackIdx];
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
    [auxDStack, auxCStack] = getStacks(spkLog, Conditions(chCond).Triggers,...
        onOffStr, timeLapse, fs, fs, spkSubs, continuousSignals);
    % Number of clusters + the piezo as the first event + the laser as the last
    % event, number of time samples in between the time window, and number of
    % total triggers.
    [Ne, Nt, NTa] = size(auxDStack);
    
    % Computing which alignment points belong to which condition.
    auxDelayFlags = false(NTa,Nccond);
    counter2 = 1;
    for ccond = consideredConditions
        auxDelayFlags(:,counter2) = ismember(Conditions(chCond).Triggers(:,1),...
            Conditions(ccond).Triggers(:,1));
        counter2 = counter2 + 1;
    end
    %% Building the population stack
    NaNew = sum(auxDelayFlags,1);
    % Removing not considered conditions
    if sum(NaNew) ~= NTa
        emtyRow = find(~sum(auxDelayFlags,2));
        auxDelayFlags(emtyRow,:) = [];
        auxDStack(:,:,emtyRow) = [];
        auxCStack(:,:,emtyRow) = [];
    end
    clInfo.id = cellfun(@(x) [sprintf('%d_',cexp), x], clInfo.id,...
        'UniformOutput', 0);
    clInfo.Properties.RowNames = clInfo.id;
    if cexp == chExp(1)
        % First assignment
        delayFlags = auxDelayFlags;
        discStack = auxDStack;
        cStack = auxCStack;
        clInfoTotal = clInfo;
        NaStack = NaNew;
    else
        % Homogenizing trial numbers
        if any(NaNew ~= NaStack)
            NaMin = min(NaStack, NaNew);
            NaMax = max(NaStack, NaNew);
            trigSubset = cell(numel(NaStack),1);
            for cc = 1:numel(NaStack)
                trigSubset{cc} = sort(randsample(NaMax(cc),NaMin(cc)));
                if NaStack(cc) == NaMin(cc)
                    tLoc = find(auxDelayFlags(:,cc));
                    tSubs = tLoc(trigSubset{cc});
                    auxDelayFlags(setdiff(tLoc,tSubs),:) = [];
                    auxDStack(:,:,setdiff(tLoc,tSubs)) = [];
                    auxCStack(:,:,setdiff(tLoc,tSubs)) = [];
                else
                    tLoc = find(delayFlags(:,cc));
                    tSubs = tLoc(trigSubset{cc});
                    delayFlags(setdiff(tLoc,tSubs),:) = [];
                    discStack(:,:,setdiff(tLoc,tSubs)) = [];
                    cStack(:,:,setdiff(tLoc,tSubs)) = [];
                    NaStack(cc) = NaNew(cc);
                end
            end
        end
        auxDStack(1,:,:) = [];
        discStack = cat(1, discStack, auxDStack);
        cStack = cat(1, cStack, auxCStack);
        popVarNames = clInfoTotal.Properties.VariableNames;
        curVarNames = clInfo.Properties.VariableNames;
        Npv = numel(popVarNames); Ncv = numel(curVarNames);
        if Npv ~= Ncv || any(~ismember(popVarNames, curVarNames))
            [logIdx, whr] = ismember(popVarNames, curVarNames);
            issueNames = find(~logIdx);
            for cvn = issueNames
                % Phy version validation
                nameFlags = contains(phyNames, popVarNames{cvn});
                if nnz(nameFlags)
                    [varRow, phyVer] = find(nameFlags);
                    [~, nPos] = ismember(curVarNames, phyNames(varRow,:));
                    nPos = logical(nPos);
                    if any(nPos)
                        fprintf(1, 'Changed %s to %s.\n',...
                            clInfo.Properties.VariableNames{nPos},...
                            popVarNames{cvn});
                        clInfo.Properties.VariableNames{nPos} =...
                            popVarNames{cvn};
                        continue
                    end
                end
                if iscell(clInfoTotal.(popVarNames{cvn}))
                    varFill = repmat({mode(string(...
                        clInfoTotal.(popVarNames{cvn})))},[size(clInfo,1),1]);
                else
                    varFill = repmat(mode(clInfoTotal.(popVarNames{cvn})),...
                        [size(clInfo,1),1]);
                end
                clInfo = addvars(clInfo, varFill, 'NewVariableNames',...
                    popVarNames{cvn});
                curVarNames = clInfo.Properties.VariableNames;
                
            end
            newVarNames = curVarNames(setdiff(1:Ncv, whr));
            if numel(newVarNames) > 1
                % If there are more than 1 variables to be included in the
                % population cluster information
                fprintf(1, '%s \n', newVarNames{:})
            elseif numel(newVarNames) == 1 &&...
                    strcmpi(newVarNames, 'amplitude')
                % If there is just one and that is 'amplitude' (very
                % specific)
                try
                    clInfo.Properties.VariableNames{setdiff(1:Ncv, whr)} =...
                        'amp';
                catch
                    clInfo = removevars(clInfo,'amp');
                    clInfo.Properties.VariableNames{setdiff(1:Ncv, whr)} =...
                        'amp';
                end
            end
        end
        clInfoTotal = cat(1, clInfoTotal, clInfo);
    end
end
gclID = clInfoTotal{logical(clInfoTotal.ActiveUnit),'id'};
Ncl = numel(gclID);
Ne = size(discStack, 1);
%% Population analysis
figureDir = fullfile(experimentDir,'PopFigures\');
if ~mkdir(figureDir)
    fprintf(1,'There was an issue with the figure folder...\n');
    fprintf(1,'Saving the figures in the data folder!\n');
    fprintf(1,'%s\n',dataDir);
    figureDir = dataDir;
end
[~,expName] = fileparts(experimentDir);
dataDir = experimentDir;
% Statistical tests
[Results, Counts] = statTests(discStack, delayFlags, timeFlags);

% Plotting statistical tests
[figs, Results] = scatterSignificance(Results, Counts,...
    consCondNames, delta_t, gclID);
arrayfun(@configureFigureToPDF, figs);
% Firing rate for all clusters, for all trials
% meanfr = cellfun(@(x) mean(x,2)/delta_t,Counts,'UniformOutput',false);
stFigBasename = fullfile(figureDir,[expName,' ']);
stFigSubfix = sprintf(' Stat RW%.1f-%.1fms SW%.1f-%.1fms',...
    responseWindow(1)*1e3, responseWindow(2)*1e3, spontaneousWindow(1)*1e3,...
    spontaneousWindow(2)*1e3);
ccn = 1;
%for cc = indCondSubs
for cc = 1:numel(figs)
    if ~ismember(cc, indCondSubs)
        altCondNames = strsplit(figs(cc).Children(2).Title.String,': ');
        altCondNames = altCondNames{2};
    else
        altCondNames = consCondNames{ccn};
        ccn = ccn + 1;
    end
    stFigName = [stFigBasename, altCondNames, stFigSubfix];
    if ~exist([stFigName,'.pdf'],'file') ||...
            ~exist([stFigName,'.emf'],'file')
        print(figs(cc),[stFigName,'.pdf'],printOpts{1}{:})
        print(figs(cc),[stFigName,'.emf'],printOpts{2})
    end
end

% Filtering for the whisker responding clusters
H = cell2mat(cellfun(@(x) x.Pvalues,...
    arrayfun(@(x) x.Activity, Results(indCondSubs), 'UniformOutput', 0),...
    'UniformOutput', 0)) < 0.05;
Htc = sum(H,2);
% Those clusters responding more than 80% of all whisker stimulating
% conditions
CtrlCond = contains(consCondNames,'control','IgnoreCase',true);
wruIdx = any(H(:,CtrlCond),2);
Nwru = nnz(wruIdx);

fprintf('%d whisker responding clusters:\n', Nwru);
fprintf('- %s\n',gclID{wruIdx})

filterIdx = true(Ne,1);
if strcmpi(filtStr, 'filtered') && nnz(wruIdx)
    filterIdx = [true; wruIdx];
end
try
clInfoTotal = addvars(clInfoTotal, false(size(clInfoTotal,1),1), 'NewVariableNames', 'Control'); 
clInfoTotal{logical(clInfoTotal.ActiveUnit), 'Control'} = wruIdx;
catch
    fprintf('Reran, perhaps?\n')
end
%% Temporal dynamics
trigTms = cell2mat(arrayfun(@(x) x.Triggers(:,1), Conditions(chCond),...
    'UniformOutput', 0)')/fs;
timesum = @(x) squeeze(sum(x,2));
clsum = @(x) squeeze(sum(x,1));
matsum = @(x) timesum(clsum(x));
stackTx = (0:Nt-1)/fs + timeLapse(1);

[~,cnd] = find(delayFlags);
[~, tmOrdSubs] = sort(trigTms, 'ascend');
cnd = cnd(tmOrdSubs); trialAx = trigTms(tmOrdSubs);
% Modulation: distance from the y=x line.
try
    modFlags = sign(clInfoTotal{clInfoTotal.ActiveUnit &...
        clInfoTotal.Control, 'Modulation'});    
catch
    clInfoTotal = addvars(clInfoTotal, zeros(size(clInfoTotal,1),1),...
        'NewVariableNames', 'Modulation');
    clInfoTotal{logical(clInfoTotal.ActiveUnit),'Modulation'} =...
        Results(1).Activity(2).Direction;
    modFlags = sign(clInfoTotal{clInfoTotal.ActiveUnit &...
        clInfoTotal.Control, 'Modulation'}); 
end

modFlags = modFlags > 0;
modFlags(:,2) = ~modFlags;
cmap = lines(Nccond);
cmap(CtrlCond,:) = ones(1,3)*1/3; cmap(:,:,2) = ones(Nccond,3)*0.7;
tdFig = figure('Name', 'Temporal dynamics', 'Color', [1,1,1]);
plotOpts = {'LineStyle', 'none', 'Marker', '.', 'Color', 'DisplayName'};
modLabel = {'Potentiation', 'Depression'};
clrSat = 0.5;
ax = gobjects(size(modFlags,2), 1);
condLey = {'Control','After Induction'};
respLey = {'Responsive', 'Non-responsive'};
% Time steps for the 3D PSTH in ms
focusStep = 2;
focusPeriods = (-2:focusStep:8)';
focusPeriods(:,2) = focusPeriods + focusStep; focusPeriods = focusPeriods * 1e-3;
Nfs = size(focusPeriods,1);
respIdx = wruIdx;
popMeanResp = zeros(Nfs, sum(NaStack), 4);

for cmod = 1:size(modFlags,2) % Up- and down-modulation
    tcount = 1;
    ax(cmod) = subplot(size(modFlags,2),1,cmod,'Parent',tdFig);
    ax(cmod).NextPlot = 'add';
    for ccond = 1:Nccond % Cycling through the conditions (control and after-induction)
        trSubs = tcount:sum(NaStack(1:ccond));
        clMod = modFlags(:,cmod);
        for cr = 1:2 % responsive and non-responsive
            for cp = 1:Nfs % 'Micro' time windows
                focusIdx = stackTx >= focusPeriods(cp,1) &...
                    stackTx <= focusPeriods(cp,2);
                clCounts = timesum(...
                    discStack( [false; respIdx], focusIdx, delayFlags(:,ccond)));
                clMean = mean(clCounts(clMod,:),1); popMean = mean(clMean);
                fprintf(1, 'Cell type: %d\n',[cr,ccond-1]*[1,2]')
                popMeanResp(cp, trSubs, [cr,cmod-1]*[1,2]') = clMean;
            end
            dispName = [condLey{ccond}, ' ', respLey{cr}, ' (',...
                num2str(popMean/diff(focusPeriod)), ' Hz)'];
            plot(ax(cmod), minutes(seconds(trialAx(trSubs))),...
                clMean/diff(focusPeriod), plotOpts{1:5}, cmap(ccond,:,cr),...
                plotOpts{6}, dispName);
            respIdx = ~respIdx;
            clMod = true(sum(respIdx),1);
        end
        tcount = 1 + sum(NaStack(1:ccond));
    end
    cmap(ccond,:,1) = brighten(cmap(ccond,:,1), clrSat);
    clrSat = -clrSat;
end

%% Firing rate distribution
%{
Ngr = pi;
while ~all(~mod(NaStack, Ngr))
    if Ngr ~= pi
        fprintf(1, '%.2f is not a divisor of %s\b!\nTry again!\n',Ngr,...
            sprintf('%d ',NaStack))
    end
    trialGroupChar = inputdlg('Trial grouping:','Trial grouping',[1,16],...
        {num2str(Ngr)});
    if str2double(trialGroupChar)
        Ngr = str2double(trialGroupChar);
    else
        disp('Cancelling the execution of the pipeline\n')
        return;
    end
end
%}

%% Getting the relative spike times for the whisker responsive units (wru)
% For each condition, the first spike of each wru will be used to compute
% the standard deviation of it.
configStructure = struct('Experiment', fullfile(dataDir,expName),...
    'Viewing_window_s', timeLapse, 'Response_window_s', responseWindow,...
    'BinSize_s', binSz, 'Trigger', struct('Name', condNames{chCond},...
    'Edge',onOffStr), 'ConsideredConditions',{consCondNames});
cellLogicalIndexing = @(x,idx) x(idx);
isWithinResponsiveWindow =...
    @(x) x > responseWindow(1) & x < responseWindow(2);

firstSpike = zeros(Nwru,Nccond);
M = 16;
binAx = responseWindow(1):binSz:responseWindow(2);
condHist = zeros(size(binAx,2)-1, Nccond);
firstOrdStats = zeros(2,Nccond);
condParams = zeros(M,3,Nccond);
txpdf = responseWindow(1):1/fs:responseWindow(2);
condPDF = zeros(numel(txpdf),Nccond);
csvSubfx = sprintf(' VW%.1f-%.1f ms.csv', timeLapse(1)*1e3, timeLapse(2)*1e3);
existFlag = false;
condRelativeSpkTms = cell(Nccond,1);
relativeSpkTmsStruct = struct('name',{},'SpikeTimes',{});
csvDir = fullfile(dataDir, 'SpikeTimes');
if ~exist(csvDir,'dir')
    iOk = mkdir(csvDir);
    if ~iOk
        fprintf(1, 'The ''Figure'' folder was not created!\n')
        fprintf(1, 'Please verify your writing permissions!\n')
    end
end
for ccond = 1:size(delayFlags,2)
    csvFileName = fullfile(csvDir,[expName,' ',consCondNames{ccond}, csvSubfx]);
    relativeSpikeTimes = getRasterFromStack(discStack,~delayFlags(:,ccond),...
        filterIdx(3:end), timeLapse, fs, true, false);
    relativeSpikeTimes(:,~delayFlags(:,ccond)) = [];
    relativeSpikeTimes(~filterIdx(2),:) = [];
    condRelativeSpkTms{ccond} = relativeSpikeTimes;
    %     respIdx = cellfun(isWithinResponsiveWindow, relativeSpikeTimes,...
    %         'UniformOutput',false);
    clSpkTms = cell(size(relativeSpikeTimes,1),1);
    if exist(csvFileName, 'file') && ccond == 1
        existFlag = true;
        ansOW = questdlg(['The exported .csv files exist! ',...
            'Would you like to overwrite them?'],'Overwrite?','Yes','No','No');
        if strcmp(ansOW,'Yes')
            existFlag = false;
            fprintf(1,'Overwriting... ');
        end
    end
    fID = 1;
    if ~existFlag
        fID = fopen(csvFileName,'w');
        fprintf(fID,'%s, %s\n','Cluster ID','Relative spike times [ms]');
    end
    rsclSub = find(filterIdx(2:end))-1;
    for cr = 1:size(relativeSpikeTimes, 1)
        clSpkTms(cr) = {sort(cell2mat(relativeSpikeTimes(cr,:)))};
        if fID > 2
            fprintf(fID,'%s,',gclID{rsclSub(cr)});
            fprintf(fID,'%f,',clSpkTms{cr});fprintf(fID,'\n');
        end
    end
    if fID > 2
        fclose(fID);
    end
    relativeSpkTmsStruct(ccond).name = consCondNames{ccond};
    relativeSpkTmsStruct(ccond).SpikeTimes = condRelativeSpkTms{ccond};
    %{
    spikeTimesINRespWin = cellfun(cellLogicalIndexing,...
        relativeSpikeTimes, respIdx, 'UniformOutput',false);
    allSpikeTimes = cell2mat(spikeTimesINRespWin(:)');
    condParams(:,:,ccond) = emforgmm(allSpikeTimes, M, 1e-6, 0);
    condPDF(:,ccond) = genP_x(condParams(:,:,ccond), txpdf);
    firstOrdStats(:,ccond) = [mean(allSpikeTimes), std(allSpikeTimes)];
    hfig = figure('Visible', 'off'); h = histogram(allSpikeTimes, binAx,...
        'Normalization', 'probability');
    condHist(:,ccond) = h.Values;
    close(hfig)
    for ccl = 1:Nwru
        frstSpikeFlag = ~cellfun(@isempty,spikeTimesINRespWin(ccl,:));
        firstSpike(ccl,ccond) = std(...
            cell2mat(spikeTimesINRespWin(ccl,frstSpikeFlag)));
    end
    %}
end
save(fullfile(csvDir,[expName,'_exportSpkTms.mat']),...
    'relativeSpkTmsStruct','configStructure')

%% ISI cumulative fraction distribution
binFig = figure('Visible','off');
areaFig = figure('Visible','on');
areaAx = axes('Parent', areaFig);
lax = log(1/fs):0.2:log(0.03);
% cmap = [0,0,102;... azul marino
%     153, 204, 255]/255; % azul cielo
cmap = lines(2);
areaOpts = {'EdgeColor', 'none', 'FaceAlpha', 0.7, 'FaceColor'};
for ccond = 1:Nccond
    ISI = cellfun(@(x) diff(x(x >= 0 & x <= 0.03)), ...
        relativeSpkTmsStruct(ccond).SpikeTimes(modFlags(:,1),:), 'UniformOutput', 0);
    ISI_merge = [ISI{:}];
    lISI = log(ISI_merge);
    hisi = histogram(lISI, lax, 'Parent', binFig, 'DisplayStyle', 'stairs');
    area(areaAx, lax(1:end-1), hisi.Values, areaOpts{:}, cmap(ccond,:));
    
    if ccond == 1
        hold(binFig.Children, 'on'); hold(areaFig.Children, 'on');
    end
end
hisi = hisi.Parent.Children;

%% Ordering PSTH
orderedStr = 'ID ordered';
dans = questdlg('Do you want to order the PSTH other than by IDs?',...
    'Order', 'Yes', 'No', 'No');
ordSubs = 1:nnz(filterIdx(2:Ncl+1));
pclID = gclID(filterIdx(2:Ncl+1));
if strcmp(dans, 'Yes')
    if ~exist('clInfoTotal','var')
        clInfoTotal = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
    end
    % varClass = varfun(@class,clInfo,'OutputFormat','cell');
    [ordSel, iOk] = listdlg('ListString', clInfoTotal.Properties.VariableNames,...
        'SelectionMode', 'multiple');
    orderedStr = [];
    ordVar = clInfoTotal.Properties.VariableNames(ordSel);
    for cvar = 1:numel(ordVar)
        orderedStr = [orderedStr, sprintf('%s ',ordVar{cvar})]; %#ok<AGROW>
    end
    orderedStr = [orderedStr, 'ordered'];
    
    if ~strcmp(ordVar,'id')
        [~,ordSubs] = sortrows(clInfoTotal(pclID,:),ordVar);
    end
end
%% Plot PSTH
goodsIdx = logical(clInfoTotal.ActiveUnit);
csNames = fieldnames(Triggers);
Nbn = diff(timeLapse)/binSz;
if (Nbn - round(Nbn)) ~= 0
    Nbn = ceil(Nbn);
end
PSTH = zeros(nnz(filterIdx) - 1, Nbn, Nccond);
for ccond = 1:Nccond
    figFileName = sprintf('%s %s VW%.1f-%.1f ms B%.1f ms RW%.1f-%.1f ms SW%.1f-%.1f ms %sset %s (%s)',...
        expName, Conditions(consideredConditions(ccond)).name, timeLapse*1e3,...
        binSz*1e3, responseWindow*1e3, spontaneousWindow*1e3, onOffStr,...
        orderedStr, filtStr);
    [PSTH(:,:,ccond), trig, sweeps] = getPSTH(discStack(filterIdx,:,:),...
        timeLapse, ~delayFlags(:,ccond), binSz, fs);
    stims = mean(auxCStack(:,:,delayFlags(:,ccond)),3);
    stims = stims - median(stims,2);
    for cs = 1:size(stims,1)
        if abs(log10(var(stims(cs,:),[],2))) < 13
            [m,b] = lineariz(stims(cs,:),1,0);
            stims(cs,:) = m*stims(cs,:) + b;
        else
            stims(cs,:) = zeros(1,Nt);
        end
    end
    figs = plotClusterReactivity(PSTH(ordSubs,:,ccond), trig, sweeps,...
        timeLapse, binSz, [{Conditions(consideredConditions(ccond)).name};...
        pclID(ordSubs)], strrep(expName,'_','\_'));
    configureFigureToPDF(figs);
    figs.Children(end).YLabel.String = [figs.Children(end).YLabel.String,...
        sprintf('^{%s}',orderedStr)];
    if ~exist([figFileName,'.pdf'], 'file') || ~exist([figFileName,'.emf'], 'file')
        print(figs, fullfile(figureDir,[figFileName, '.pdf']),'-dpdf','-fillpage')
        print(figs, fullfile(figureDir,[figFileName, '.emf']),'-dmeta')
    end
end

%% Comparing the PSTHs for all conditions
psthFig = figure('Color', [1,1,1], 'Name', 'Condition PSTH comparison',...
    'Visible', 'off');
txpsth = (0:Nbn-1)*binSz + timeLapse(1);
focusWindow = [-2, 30]*1e-3;
focusIdx = txpsth >= focusWindow(1) & txpsth <= focusWindow(2);
txfocus = txpsth(focusIdx);
PSTH_raw = squeeze(sum(PSTH(modFlags(:,2),:,:),1));
PSTH_trial = PSTH_raw./NaStack;
PSTH_prob = PSTH_raw./sum(PSTH_raw,1);
PSTH_all = cat(3, PSTH_raw, PSTH_trial, PSTH_prob);
axp = gobjects(3,1);
subpltsTitles = {'PSTH for all conditions (trial-normalized)',...
    'Cumulative density function', 'Cumulative sum for normalized spikes'};
yaxsLbls = {'Spikes / Trials', 'Spike probability', 'Spike number'};
for cax = 1:3
    axp(cax) = subplot(4, 1, cax, 'Parent', psthFig, 'NextPlot', 'add');
    for ccond = 1:Nccond
        plotOpts = {'DisplayName', consCondNames{ccond}};
        switch cax
            case 1 % Plot the trial-normalized PSTH; page #2
                plot(axp(cax), txfocus, PSTH_all(focusIdx,ccond,2),...
                    plotOpts{:})
            case 2 % Plot the cumulative probability function; page 3
                PSTH_aux = PSTH_all(focusIdx, ccond, 3);
                PSTH_aux = PSTH_aux / sum(PSTH_aux);
                plot(axp(cax), txfocus, cumsum(PSTH_aux), plotOpts{:})
            case 3 % Plot the cumulative sum for the mean spikes; page 2
                plot(axp(cax), txfocus, ...
                    cumsum(PSTH_all(focusIdx, ccond, 2)), plotOpts{:})
        end
    end
    legend(axp(cax), 'show', 'Location', 'best')
    title(axp(cax), subpltsTitles{cax})
    ylabel(axp(cax), yaxsLbls{cax})
    if cax < 3
        xticks(axp(cax), '')
    else
        xlabel(axp(cax), sprintf('Time_{%.2f ms} [s]', binSz*1e3))
    end
end
axp(4) = subplot(4,1,4, 'Parent', psthFig);
PSTH_diff = (PSTH_trial(:,2) - PSTH_trial(:,1))./PSTH_trial(:,1);
bar(axp(4), txfocus(PSTH_diff(focusIdx) > 0), PSTH_diff(PSTH_diff > 0 & focusIdx'),...
    'FaceColor', [51, 204, 51]/255, 'DisplayName', 'Potentiation'); hold on
bar(axp(4), txfocus(PSTH_diff(focusIdx) <= 0), PSTH_diff(PSTH_diff <= 0 & focusIdx'),...
    'FaceColor', [204, 51, 0]/255, 'DisplayName', 'Depression');
axp(4).Box = 'off';
linkaxes(axp, 'x')
legend show
psthFig.Visible = 'on';
psthFig = configureFigureToPDF(psthFig);
psthFigFileName = sprintf('%s PSTH All conditions RW%.2f-%.2f ms',expName,...
    responseWindow*1e3);
psthFigFileName = fullfile(figureDir,psthFigFileName);
if ~exist([psthFigFileName,'.pdf'], 'file') 
    print(psthFig, [psthFigFileName, '.pdf'],'-dpdf','-fillpage')
end
if ~exist([psthFigFileName,'.emf'], 'file')
    print(psthFig, [psthFigFileName, '.emf'],'-dmeta')
end


%%