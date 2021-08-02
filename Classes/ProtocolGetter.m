classdef ProtocolGetter < handle
    %PROTOCOLGETTER tries to recognize and extract the experimental
    %protocol from either single or merged recordings
    
    properties (SetAccess = 'private')
        dataDir char = '';
        ismerged (1,1) logical = false;
        fileOrder string = "";
        condSigFiles string = "";
        fs (1,1) double = 3e4;
        Triggers struct;
        Edges (1,2) struct;
        Conditions struct;
        isSaved logical = false;
        BinFile string = "";
    end
    
    properties (SetAccess = 'private',GetAccess = 'protected')
        awaken logical = false;
    end
    
    methods
        function obj = ProtocolGetter(directory)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            foStr = '_fileOrder';
            getBaseNames = @(x) arrayfun(@(y) fileparts(y.name), x,...
                'UniformOutput', 0);
            
            % Good folder?
            if ~exist(directory,'dir')
                fprintf('Invalid folder! No object created!\n')
                return
            end
            obj.dataDir = directory;
            fsFile = dir(fullfile(directory,'*_sampling_frequency.mat'));
            load(fullfile(directory,fsFile(1).name),'fs')
            obj.fs = fs;
            [~,binBaseNames] = getBaseNames(dir(fullfile(directory, '*.bin')));
            smrxFiles = dir(fullfile(directory, '*.smrx'));
            subfolderFlag = false;
            if isempty(smrxFiles)
                subfolderFlag = true;
                smrxFiles = dir(fullfile(directory,'*\*.smrx'));
            elseif isempty(smrxFiles)
                fprintf(1,'The given directory contains no .smrx files. Please ')
                fprintf(1,'try again with another folder which do contain .smrx files.\n')
                return
            end
            [~,smrxBaseNames] = getBaseNames(smrxFiles);
            merFiles = dir(fullfile(directory, ['*', foStr, '.txt']));
            Nb = numel(binBaseNames); Ns = numel(smrxBaseNames);
            % Checking if the number of binnary files equals the number of
            % CED recorded files and if their names match. sitFlags has 4
            % likely cases encoded from.
            sitFlags = [Nb == Ns, isempty(merFiles),...
                any(ismember(binBaseNames, smrxBaseNames))];
            selNum = bin2num(sitFlags,0);
            if any(selNum == setdiff(0:5,3))
                if ~sitFlags(2)
                    obj.ismerged = true;
                    foID = fopen(fullfile(merFiles(1).folder,...
                        merFiles(1).name),'r');
                    smrxFlags = zeros(Ns,1,'single'); counter = 1;
                    smrxBaseNamesStr = string(smrxBaseNames);
                    while ~feof(foID) && counter <= Ns
                        [~,cfName,~] = fileparts(string(fgetl(foID)));
                        smrxFlags = smrxFlags + counter * strcmpi(...
                            smrxBaseNamesStr, cfName);
                        counter = counter + 1;
                    end
                    smrxBaseNames(smrxFlags == 0) = [];
                    smrxFlags(smrxFlags == 0) = [];
                    if feof(foID)
                        binFile = cfName + ".bin";
                    else
                        binFile = string(fgetl(foID));
                    end
                    fclose(foID);
                    if subfolderFlag
                        [~,subFolders] = arrayfun(@(x) fileparts(x.folder),...
                                    smrxFiles, 'UniformOutput', 0);
                        smrxBaseNames(smrxFlags) = arrayfun(@(x,y) ...
                            fullfile(x,y), subFolders, smrxBaseNames);
                    else
                        smrxBaseNames(smrxFlags) = smrxBaseNames;
                    end
                else
                    % User interaction for older merges
                    mergeQuest =...
                        sprintf(...
                        ['Did you merge %s.bin from more than one ''smrx''',...
                        ' file?'],binBaseNames{1});
                    mm = questdlg(mergeQuest,...
                        'Merged?','Yes','No','Yes');
                    % Auxiliary variable for awaken experiments
                    binRecordFile = 1;
                    if strcmpi(mm,'yes')
                        % Which smrx files were used
                        obj.ismerged = true;
                        [smrxFlags, iOk] = listdlg(...
                            'ListString',smrxBaseNames,...
                            'SelectionMode','multiple',...
                            'Name','SMRX select',...
                            'PromptString','Out of which smrx files?');
                        if iOk
                            % Correct order?
                            smrxBaseNames = smrxBaseNames(smrxFlags);
                            oo = questdlg('Are they in order?',...
                                'Order','Yes','No','Yes');
                            if strcmpi(oo,'no')
                                % Reorder the file names
                                orderAns = inputdlg(smrxBaseNames,...
                                    'File order',[1, 60],string(smrxFlags));
                                nFileOrder = str2double(orderAns);
                                nSmrxBaseNames = smrxBaseNames;
                                if ~isempty(orderAns) && ~issorted(nFileOrder)
                                    fprintf(1,'Changing file order...\n')
                                    nSmrxBaseNames(nFileOrder) = smrxBaseNames;
                                    smrxBaseNames = nSmrxBaseNames;
                                    smrxFlags = nFileOrder;
                                else
                                    fprintf(1, 'Order unaltered.\n')
                                end
                                
                            end
                            % Create a merge file (_fileOrder)
                            binSl = 1;
                            if Nb > 1
                                [binSl, iOk] =...
                                    listdlg('ListString',binBaseNames,...
                                    'SelectionMode','single',...
                                    'PromptString','For which binary file');
                                if ~iOk
                                    fprintf(1, 'Error: None binary file selected!\n');
                                    fprintf(1, 'No object created.\n')
                                    return
                                end
                            end
                             subFolders = '';
                            if subfolderFlag
                                [~,subFolders] = arrayfun(...
                                    @(x) fileparts(x.folder),...
                                    smrxFiles,...
                                    'UniformOutput', 0);
                                subFolders(nFileOrder) = subFolders;
                                subFolders = string(subFolders);
                                smrxBaseNames = arrayfun(@(x,y) fullfile(...
                                    x,y), subFolders, smrxBaseNames);
                            end
                            foID = fopen(fullfile(obj.dataDir,...
                                [binBaseNames{binSl}, '_fileOrder.txt']),'w');
                            for csf = 1:size(smrxBaseNames,1)
                                fprintf(foID, '%s\n',...
                                    [smrxBaseNames{csf}, '.smrx']);
                            end
                            fprintf(foID, '%s', [binBaseNames{binSl},...
                                '.bin']);
                            fclose(foID);
                        else
                            fprintf(1, 'Cancelled by the user.\n');
                            fprintf(1, 'No object created.\n');
                            return
                        end
                    else
                        % Binary file builed from one smrx
                        smrxFlags = 1;
                        if Ns > 1 && sitFlag(3)
                            % If there are more than one but one matching
                            % name
                            smrxUniqueFlag = ismember(smrxBaseNames,...
                                binBaseNames(1));
                            if smrxUniqueFlag
                                smrxBaseNames = smrxBaseNames(smrxUniqueFlag);
                            else
                                % No matching name
                                [usedSxFile, iOk] =...
                                    listdlg('ListString',smrxBaseNames,...
                                    'Name','SMRX file selection',...
                                    'PromptString','Select the SMRX file',...
                                    'SelectionMode','single');
                                if iOk
                                    smrxBaseNames = smrxBaseNames(usedSxFile);
                                else
                                    % File not in this folder?
                                    fprintf(1,['Perhaps the file is not',...
                                        'located in this folder?\n',...
                                        'No object created.\n'])
                                    return
                                end
                            end
                        else
                            % Didn't record with Spike2, but with other
                            % program directly to bin (most likely)
                            awakeAns = questdlg(['Did you convert an SMRX',...
                                ' file to binary?'],'SMRX used?',...
                                'Yes','No','No');
                            if strcmpi(awakeAns,'Yes')
                                fprintf(1,['Couldn''t figure out which ',...
                                    'file is the original nor the trigger!',...
                                    ' Cannot continue...\nNo object created!\n'])
                            end
                            obj.awaken = true;
                            if Nb > 1
                                [binRecordFile, iOk] = listdlg('ListString',...
                                    binBaseNames, 'SelectionMode', 'multiple');
                                if ~iOk
                                    fprintf(1, 'Cancelling...\nNo object created!\n')
                                end
                            end
                        end
                    end
                    binFile = binBaseNames{binRecordFile} + ".bin";
                end
                if obj.awaken
                    obj.fileOrder = string(binBaseNames{binRecordFile}) + ".bin";
                else
                    obj.fileOrder = string(smrxBaseNames) + ".smrx";
                end
            else
                % Likely single file transformation
                binFile = binBaseNames{1} + ".bin";
                obj.fileOrder = string(smrxBaseNames) + ".smrx";
            end
            obj.BinFile = binFile;
        end
        
        function obj = getConditionSignals(obj)
            %GETCONDITIONSIGNALS extracts the condition signals from the
            %SMRX files or from the TriggerSignals*.bin file.
            if obj.awaken
                fprintf(1,'Reading BIN file:\n');
                trigFile = dir(fullfile(obj.dataDir,'TriggerSignals*.bin'));
                if ~isempty(trigFile)
                    fID = fopen(fullfile(trigFile.folder,trigFile.name),'r');
                    trig = fread(fID, Inf, 'uint16'); [~] = fclose(fID);
                    trig = reshape(trig, 2, []); trig = trig - 2^15;
                    trig = int16(trig); trig = trig - median(trig, 2);
                    trig(1,:) = -trig(1,:);
                    baseName = strsplit(obj.BinFile,'.bin'); baseName = baseName(1);
                    condSigFile = fullfile(obj.dataDir,string(baseName)+"_AwakeCondSig.mat");
                    save(condSigFile, 'trig'); obj.condSigFiles = condSigFile;
                end
            else
                fprintf(1,'Reading SMRX files:\n')
                Ns = size(obj.fileOrder,1);
                for csf = 1:Ns
                    fprintf(1,'%s\n',obj.fileOrder(csf));
                    fID = fopen(fullfile(obj.dataDir, obj.fileOrder(csf)),'r');
                    getConditionSignalsBF(fID);
                end
                obj.condSigFiles = extractBefore(obj.fileOrder,".smrx") +...
                    "_CondSig.mat";
            end
        end
        
        function obj = getSignalEdges(obj)
            %GETSIGNALEDGES looks into the _condSig files and searches for
            %whisker, laser, and lfp channels for then concatenate the edge
            %samples.
            mStimStruct = struct();
            for csf = 1:size(obj.condSigFiles,1)
                % Load the channel variables
                stimSig = load(...
                    fullfile(obj.dataDir,obj.condSigFiles(csf)));
                fields = fieldnames(stimSig);
                % Identify who is who
                [idMat, titles, headers] =...
                    ProtocolGetter.searchIDFromSignals(stimSig);
                headers = string(headers);
                % Variable assignment
                stSgStruct = ProtocolGetter.assign2StimulationSignals(...
                    stimSig, idMat, titles, fields);
                wHead = stimSig.(headers(idMat(:,1)));
                try
                    eHead = stimSig.(headers(idMat(:,3)));
                    stSgStruct = ProtocolGetter.correctLFPLength(...
                        stSgStruct, wHead, eHead);
                catch
                    fprintf(1, 'No LFP in the current experiment!\n')
                end
                mStimStruct = catStruct(mStimStruct, stSgStruct);
            end
            obj.Triggers = mStimStruct;
            % Extract signal edges
            [wSubs, lSubs] = ProtocolGetter.extractSignalEdges(mStimStruct,...
                obj.fs);
            subs = {wSubs,lSubs};
            for cstr = 1:2
                obj.Edges(cstr).Name = titles{idMat(:,cstr)};
                obj.Edges(cstr).Subs = subs{cstr};
            end
            function mStimStruct = catStruct(mStimStruct, stSgStruct)
                fns = string(fieldnames(mStimStruct));
                if isempty(fns)
                    mStimStruct = stSgStruct;
                else
                    for cfn = 1:length(fns)
                        if ~isempty(stSgStruct.(fns(cfn)))
                            mStimStruct.(fns(cfn)) = cat(...
                                2*isrow(stSgStruct.(fns(cfn))) +...
                                1*iscolumn(stSgStruct.(fns(cfn))),...
                                mStimStruct.(fns(cfn)),stSgStruct.(fns(cfn)));
                        end
                    end
                end
            end
        end
        
        function obj = getFrequencyEdges(obj)
            %GETFREQUENCYEDGES looks for train of stimulus in both signals
            %and saves the position and frequency of the train.
            % Whisker, piezo, mechanical
            for cst = 1:size(obj.Edges, 2)
                coupFreq = 0;
                % Frequencies found (wFreq), First pulse of the trains
                % (wFlags), and frequency by pulse (wFreqs)
                [wFreq, wFlags, wFreqs] = ProtocolGetter.extractFrequencyTrains(...
                    obj.Edges(cst).Subs, obj.fs);
                wTrainBodyFlags = ismembertol([0;wFreqs], wFreq, 0.01);
                obj.Edges(cst).Subs = obj.Edges(cst).Subs(~wTrainBodyFlags,:);
                wFlags(wTrainBodyFlags) = [];
                wFreqs(wTrainBodyFlags(1:end-1)) = [];
                if ~isempty(wFreq)
                    coupFreq = round(wFlags .* wFreqs, 1);
                end
                obj.Edges(cst).Frequency = coupFreq;
                obj.Edges(cst).FreqValues = wFreq;
                fprintf(1, 'Found %d frequencies for %s (', numel(wFreq),...
                    obj.Edges(cst).Name)
                cf = 1;
                while cf <= numel(wFreq)
                    if cf < numel(wFreq)
                        fprintf(1, '%.1f, ', wFreq(cf))
                    else
                        fprintf(1, '%.1f', wFreq(end))
                    end
                    cf = cf + 1;
                end
                fprintf(1, ' Hz)\n')
            end
        end
        
        function obj = pairStimulus(obj)
            %PAIRSTIMULUS looks at the temporal relationships between both
            %stimulus. In other words, the delay between each other.
            % Regardless of the pairing or condition, all stimuli from one
            % signal are going to be grouped in the first conditions and
            % labeled as '<signal_name>All', replacing <signal_name> by the
            % actual name of your signal.
            fetchSubs = @(idx, bf, subOrd, subs) ...
                subs(sort(subOrd(bf(:,idx))),:);
            wSub = obj.Edges(1).Subs;lSub = obj.Edges(2).Subs;
            wFreqs = obj.Edges(1).Frequency;lFreqs = obj.Edges(2).Frequency;
            wFreq = obj.Edges(1).FreqValues;lFreq = obj.Edges(2).FreqValues;
            validateFreq = @(idx, bf, subOrd, freqs)...
                nnz(fetchSubs(idx, bf, subOrd, freqs));
            removeZeros = @(x) x(x~=0);
            for ccon = 1:size(obj.Edges,2)  %#ok<*FXUP>
                obj.Conditions(ccon).name = [obj.Edges(ccon).Name, 'All'];
                obj.Conditions(ccon).Triggers = obj.Edges(ccon).Subs;
            end
            % Finding the time difference between every pulse in two
            % signals
            mxPulses = min(size(lSub,1),size(wSub,1));
            dm = distmatrix(lSub(:,1)/obj.fs,wSub(:,1)/obj.fs);
            [strDelay, whr] = sort(dm(:),'ascend');
            [lSubOrd, wSubOrd] = ind2sub(size(dm),whr(1:mxPulses));
            timeDelay = strDelay(1:mxPulses);
            % The delay will usually be milliseconds long, so a logarithmic
            % scale will be useful.
            delays = 10.^uniquetol(log10(timeDelay),0.01/log10(max(abs(timeDelay))));
            % Removing delays that are greater than 1 second.
            delays(delays > 1) = [];
            if std(delays.*1e3) < 1
                % Validation for similarity between the delays. If the
                % standard deviation of the delays is smaller than 1 ms,
                % then it is very likely that there are no protocolled
                % delays.
                delays = mean(delays);
            end
            % Assigning the subscripts to the Condition structure.
            % Total number of delays
            Ndel = numel(delays);
            fprintf(1,'%d Delays found:', Ndel)
            % Logical matrix indicating membership of the subscripts to one
            % or the other delays.
            lsDel = false(length(timeDelay),Ndel);
            Ncond = numel(obj.Conditions);
            tol = (1e-4)/max(timeDelay);
            for cdl = 1:Ndel
                % Starting from the last condition on
                fprintf(1,' %.1f',delays(cdl)*1e3)
                % Assign the boolean membership
                % lsDel(:,cdl) = ismembertol(log10(timeDelay),log10(delays(cdl)),...
                %     abs(0.01/log10(max(delays))));
                lsDel(:, cdl) = ismembertol(timeDelay, delays(cdl), tol);
                % Create the name of the condition
                obj.Conditions(Ncond + cdl).name = sprintf('Delay %0.3f s',...
                    delays(cdl));
                % Verify if the selected subscripts are pulse trains 
                if (all(wFreqs) || all(lFreqs)) &&...
                        ((validateFreq(cdl, lsDel, lSubOrd, lFreqs) ||...
                        validateFreq(cdl, lsDel, wSubOrd, wFreqs)))
                    delLFreqs = removeZeros(unique(...
                        fetchSubs(cdl, lsDel, lSubOrd, lFreqs)));
                    delWFreqs = removeZeros(unique(...
                        fetchSubs(cdl, lsDel, wSubOrd, wFreqs)));
                    for cdf = 1:length(delLFreqs)
                        obj.Conditions(Ncond + cdl).name =...
                            [obj.Conditions(Ncond + cdl).name,...
                            sprintf(' + L%.1f',delLFreqs(cdf))];
                    end
                    allDelays = sort(wSubOrd(lsDel(:,cdl)));
                    for cdf = 0:length(delWFreqs)-1
                        % Adding different frequencies for the given delay
                        freqSubIdx = ismember(allDelays,...
                            find(ismember(wFreqs,delWFreqs(cdf+1))));
                        obj.Conditions(Ncond + cdl + cdf).Triggers = ...
                            wSub(allDelays(freqSubIdx),:);
                        obj.Conditions(Ncond + cdl + cdf).name =...
                            [sprintf('Delay %0.3f s', delays(cdl)),...
                            sprintf(' + W%.1f',delWFreqs(cdf+1))];
                    end
                    continue
                end
                % Use the boolean membership to find the subscripts that
                % belong to the current condition.
                obj.Conditions(Ncond + cdl).Triggers =...
                    fetchSubs(cdl, lsDel, wSubOrd, wSub);
                %  wSub(sort(wSubOrd(lsDel(:,cdl))),:); Line before
            end
            Ncond = numel(obj.Conditions);
            fprintf(1, ' ms\n')
            delFlag = any(lsDel,2);
            % Removing the used subscripts for delays (and possibly some
            % frequencies)
            wSub(fetchSubs(1,delFlag,wSubOrd,(1:size(wSub,1))'),:) = [];
            if nnz(wFreqs)
                wFreqs(fetchSubs(1,delFlag,wSubOrd,(1:size(wFreqs,1))')) = [];
            end
            lSub(fetchSubs(1,delFlag,lSubOrd,(1:size(lSub,1))'),:) = [];
            if nnz(lFreqs)
                lFreqs(fetchSubs(1,delFlag,lSubOrd,(1:size(lFreqs,1))')) = [];
            end
            % Creating the unpaired frequency conditions and removing those
            % subscripts
            [obj, wSub] = addFrequencyStimulus(obj,wSub,wFreqs,wFreq,1);
            [obj, lSub] = addFrequencyStimulus(obj,lSub,lFreqs,lFreq,2);
            % Add finally the control conditions
            obj.Conditions(Ncond + 1).name =['Control ',obj.Edges(1).Name];
            obj.Conditions(Ncond + 1).Triggers = wSub;
            obj.Conditions(Ncond + 2).name =['Control ',obj.Edges(2).Name];
            obj.Conditions(Ncond + 2).Triggers = lSub;
            function [obj, subs] =...
                    addFrequencyStimulus(obj,subs,freqs,freqVals,idx)
                [Subf, FrFlags, CondNames] = fetchFrequencySubs(subs,...
                    freqs, freqVals);
                Nnc = size(Subf,1);
                for ccon = 1:Nnc
                    obj.Conditions(Ncond + ccon).name = [obj.Edges(idx).Name,...
                        ' ', CondNames{ccon}];
                    obj.Conditions(Ncond + ccon).Triggers = Subf{ccon};
                end
                subs(any(FrFlags,2),:) = [];
                Ncond = numel(obj.Conditions);
                function [fsubs, msFlag, names] =...
                        fetchFrequencySubs(subs, freqs, freqVals)
                    Nf = numel(freqVals);
                    fsubs = cell(Nf,1);
                    names = cell(fsubs);
                    msFlag = false(size(freqs,1),Nf);
                    for cf = 1:Nf
                        msFlag(:,cf) = ismembertol(freqs,freqVals(cf),...
                            0.11, 'DataScale', 1);
                        fsubs{cf} = subs(msFlag(:,cf),:);
                        names{cf} = sprintf('%.1f Hz', freqVals(cf));
                    end
                end
            end
        end
        
        function obj = saveConditions(obj)
            if ~obj.isSaved
                binBaseName = obj.BinFile.extractBefore(".bin");
                if ~isempty(obj.Conditions) && ~isempty(obj.Triggers)
                    Conditions = obj.Conditions; %#ok<*NASGU,*PROP>
                    Triggers = obj.Triggers; fs = obj.fs;
                    condFileName = fullfile(obj.dataDir, binBaseName+...
                        "analysis.mat");
                    save(condFileName, 'Conditions', 'Triggers', 'fs')
                    obj.isSaved = true;
                end
            end
        end
        
    end
    
    methods (Static, Access = 'private')
        function [idMat, titles, headers] = searchIDFromSignals(stimSig)
            checkSignal = @(x,y) contains(x,y,'IgnoreCase',true);
            fields = fieldnames(stimSig);
            chanFlag = cellfun(@contains,fields,repmat({'chan'},numel(fields),1));
            chanSubs = find(chanFlag);
            headers = cellfun(@strrep,fields(chanFlag),...
                repmat({'chan'},numel(chanSubs),1),...
                repmat({'head'},numel(chanSubs),1),'UniformOutput',false);
            titles = cell(numel(headers),1);
            whiskFlag = false(numel(titles),1);
            laserFlag = whiskFlag;
            lfpFlag = whiskFlag;
            for chead = 1:numel(headers)
                titles{chead} = stimSig.(headers{chead}).title;
                whiskFlag(chead) = checkSignal(titles{chead},'piezo') |...
                    checkSignal(titles{chead},'puff') |...
                    checkSignal(titles{chead},'mech');
                laserFlag(chead) = checkSignal(titles{chead},'laser');
                lfpFlag(chead) = checkSignal(titles{chead},'lfp');
            end
            idMat = [whiskFlag, laserFlag, lfpFlag];
        end
        
        function stSgStruct =...
                assign2StimulationSignals(stimSig, idMat, titles, fields)
            whiskSubs = 1:numel(titles);
            laserSubs = whiskSubs;
            lfpSubs = whiskSubs;
            whiskFlag = idMat(:,1); laserFlag = idMat(:,2);
            lfpFlag = idMat(:,3);
            chanSubs = find(contains(fields,'chan'));
            if any(whiskFlag)
                whiskSubs = find(whiskFlag);
            end
            while sum(whiskFlag) ~= 1
                wSub = listdlg('ListString',titles(whiskSubs),...
                    'PromptString','Select the mechanical TTL',...
                    'SelectionMode','single');
                if ~isempty(wSub)
                    whiskFlag = false(size(whiskFlag));
                    whiskFlag(wSub) = true;
                else
                    fprintf(1,'Please select one of the displayed signals!\n')
                end
            end
            whisk = stimSig.(fields{chanSubs(whiskFlag)});
            if any(laserFlag)
                laserSubs = find(laserFlag);
            end
            while sum(laserFlag) ~= 1
                lSub = listdlg('ListString',titles(laserSubs),...
                    'PromptString','Select the laser TTL',...
                    'SelectionMode','single');
                if ~isempty(lSub)
                    laserFlag = false(size(laserFlag));
                    laserFlag(lSub) = true;
                else
                    fprintf(1,'Please select one of the displayed signals!\n')
                end
            end
            laser = stimSig.(fields{chanSubs(laserFlag)});
            
            if any(lfpFlag)
                lfpSubs = find(lfpFlag);
            end
            iOk = true;
            while sum(lfpFlag) ~= 1 && iOk
                [lSub, iOk] = listdlg('ListString',titles(lfpSubs),...
                    'PromptString','Select the lfp signal',...
                    'SelectionMode','single','CancelString', 'none');
                if ~isempty(lSub)
                    lfpFlag = false(size(lfpFlag));
                    lfpFlag(lSub) = true;
                end
            end
            if iOk
                lfp = stimSig.(fields{chanSubs(lfpFlag)});
            else
                lfp = 0;
            end
            stSgStruct = struct('Whisker',whisk,'Laser',laser,'LFP',lfp);
        end
        
        function [stSgStruct, Ns] = correctLFPLength(stSgStruct, wHead,...
                eHead)
            expDomain = round([min(wHead.start, eHead.start),...
                max(wHead.stop, eHead.stop)] *...
                max(wHead.SamplingFrequency, eHead.SamplingFrequency));
            intanFs = wHead.SamplingFrequency;
            intanDomain = round([wHead.start, wHead.stop]*wHead.SamplingFrequency);
            cedFs = eHead.SamplingFrequency;
            cedDomain = round([eHead.start, eHead.stop] * cedFs);
            tx = (0:expDomain(2)-1)/intanFs;
            lfp = stSgStruct.LFP;
            if cedFs ~= intanFs
                lfp = resample(double(lfp), intanFs/cedFs, 1);
                cedDomain = cedDomain .* [1, intanFs/cedFs];
            end
            
            % Length LFP  correction
            init = max(intanDomain(1), cedDomain(1));
            if cedDomain(2) < intanDomain(2)
                fin = cedDomain(2);
            else
                fin = intanDomain(2);
            end
            stSgStruct.LFP = lfp(init:fin);
            
            % Intan signals corrections
            init = 1;
            fin = size(stSgStruct.LFP);
            flds = fieldnames(stSgStruct);
            lfpFlag = contains(flds, 'lfp', 'IgnoreCase', 1);
            flds(lfpFlag) = [];
            for cf = 1:numel(flds)
                vars = stSgStruct.(flds{cf});
                stSgStruct.(flds{cf}) = vars(init:fin);
            end
            Ns = fin;
        end
        
        function [wSubs, lSubs] = extractSignalEdges(Triggers, fs)
            wObj = StepWaveform(Triggers.Whisker, fs);
            lObj = StepWaveform(Triggers.Laser, fs);
            getGlitch = @(x) (diff(x,1,2) == 0);
            wSubs = wObj.subTriggers;
            wSubs(getGlitch(wSubs),:) = [];
            lSubs = lObj.subTriggers;
            lSubs(getGlitch(lSubs),:) = [];
        end
        
        function [freqCond, fstSubs, pulsFreq] =...
                extractFrequencyTrains(subs, fs)
            
            tms = subs(:,1)/fs;
            pulsFreq = 1./diff(tms);
            % Logical index pointing at the first pulse of a frequency
            % train
            fstSubs = StepWaveform.firstOfTrain(tms);
            % Inverse of the time difference between pulses (frequency)
            
            freqCond = round(uniquetol(pulsFreq, 0.1/max(pulsFreq)), 1);
            freqCond = freqCond(freqCond >= 1); % Empty for no frequency
        end
        
    end
end