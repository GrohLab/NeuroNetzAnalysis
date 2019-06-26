classdef UMSDataLoader < handle
    % DATALOADER
    properties (SetAccess = 'private')
        Data(1,:) cell = {};
        SpikeTimes
        SpikeUMSStruct;
    end
    
    properties
        SamplingFrequency (1,1) double = 2e4;
        % Detecting parameters
        DetectMethod char = 'auto';
        Thresh (1,1) double= 3.9;
        Shadow (1,1) double= 0.9;             % ms, enforced dead region after each spike
        RefractoryPeriod (1,1) double= 2.5;   % ms, refractory period (for calculation refractory period violations)
        % Alignment parameters
        WindowSize (1,1) double = 1.5;        % ms, width of a spike
        CrossTime (1,1) double = 0.6;         % ms, alignment point for peak of waveform
        MaxJitter (1,1) double = 0.6;         % ms, width of window used to detect peak after threshold crossing
        % sorting parameters
        AggCutoff (1,1) double = .05;         %  higher = less aggregation, lower = more aggregation
        KmeansClusterSize (1,1) double = 500; %  target size for miniclusters
    end % properties
    properties (Dependent)
        Ns                      % Number of samples
        Nch                     % Number of channels
        Time                    % Time axis
    end
    properties (SetAccess = 'private',GetAccess = 'private')
        % Default value is for the Poly design unaccesible to the user:
        PolyChanOrder (1,:) int16 =...
            [8, 9, 7, 10, 4, 13, 5, 12, 2, 15, 1, 16, 6, 11, 3, 14];
        FileName char;
        ChannelsID;
    end
    methods
        function obj = UMSDataLoader(arg1,arg2)
            % arg1 could be either a filename or a numeric vector .
            % arg2 could be either the sampling frequency or the channel
            % number.
            % For UMS2k to process the spike traces, it needs a cell array
            % of Nt (number of trials) which contains a Ns x Nch (number of
            % samples x number of channels) matrix. This object arranges
            % the given data into such format from a file or from the
            % workspace.
            if exist('arg1','var') && ischar(arg1)
                if nargin == 1
                    [data, fs, chanReadOut, fileNameOut, chansID] =...
                        UMSDataLoader.LoadFromFile(arg1, []);
                elseif nargin == 2
                    [data, fs, chanReadOut, fileNameOut, chansID] =...
                        UMSDataLoader.LoadFromFile(arg1, arg2);
                end
            elseif exist('arg1','var') && isnumeric(arg1)
                if isrow(arg1)
                    arg1 = arg1';
                end
                data = {arg1};
                if exist('arg2','var') && isnumeric(arg2)
                    obj.SamplingFrequency = arg2;
                    fs = arg2;
                end
                chanReadOut = 1;
                fileNameOut = 'Data_from_workspace';
                chansID = 1;
            end
            if nargin == 0
                data = [];
            end
            if ~isempty(data)
                obj.FileName = fileNameOut;
                obj.Data = data;
                obj.changeChannelOrder(chanReadOut);
                obj.SamplingFrequency = fs;
                obj.ChannelsID = chansID;
            end
            
            disp('Constructed!')
        end % Constructor
        
        function invertSignals(obj)
            obj.Data{1} = -obj.Data{1};
        end
        
        function data = getDataMatrix(obj)
            data = obj.Data{1};
        end
        
        %% Ultra Mega Sort 2000 Pipeline
        function obj = UMS2kPipeline(obj)
            spikesLocal = ss_default_params(obj.SamplingFrequency,...
                'thresh',obj.Thresh,...
                'window_size',obj.WindowSize,'shadow',obj.Shadow,...
                'cross_time',obj.CrossTime,...
                'refractory_period',obj.RefractoryPeriod,...
                'max_jitter',obj.MaxJitter,'agg_cutoff',obj.AggCutoff,...
                'kmeans_clustersize',obj.KmeansClusterSize);
            % The data is ordered as it is given to UMS2k to process.
            % The ordering of the channels is not necessary at this point.
            % Probably a good idea is to order the channels already by the
            % beguinning.
            tempData{1} = obj.Data{1};
            spikesLocal = ss_detect(tempData, spikesLocal);
            spikesLocal = ss_align(spikesLocal);
            spikesLocal = ss_kmeans(spikesLocal);
            spikesLocal = ss_energy(spikesLocal);
            spikesLocal = ss_aggregate(spikesLocal);
            splitmerge_tool(spikesLocal)
            h = gcf;
            
            % getSpikeStructureCls is a private method which is only used
            % in this context.
            set(h,'CloseRequestFcn',{@obj.getSpikeStructureCls,h})
            
            % When the figure is closed, the spikes structure is saved in
            % the object for spike extraction from the human recognised
            % cluster.
            try
                while strcmp(h.BeingDeleted,'off')
                    waitforbuttonpress
                end
            catch
            end
            
        end
        
        function changeDataFromFile(obj,fileName,chanOrder)
            % Loads the channels with the provided numbers
            if nargin == 2
                fprintf('Importing all channels in their original order.\n')
                [data, fs, chanReadOut, fileNameOut, channIDs] =...
                    obj.LoadFromFile(fileName,[]);
            elseif nargin == 3
                [data, fs, chanReadOut, fileNameOut, channIDs] =...
                    obj.LoadFromFile(fileName,chanOrder); 
            end
            if ~isempty(data)
                obj.Data = data;
                obj.SamplingFrequency = fs;
                obj.changeChannelOrder(chanReadOut);
                obj.FileName = fileNameOut;
                obj.ChannelsID = channIDs;
                obj.SpikeTimes = [];
                obj.SpikeUMSStruct = [];
            end
        end
        
        %% SET AND GET SpikeUMSStruct
        function changeUMSStructure(obj,structIn)
            obj.SpikeUMSStruct = structIn;
            obj.SamplingFrequency = structIn.params.Fs;
        end
        
        function set.SpikeUMSStruct(obj,structIn)
            if isstruct(structIn) && isfield(structIn,'params')
                obj.SpikeUMSStruct = structIn;
            else
                fprintf('Deleting the UMS spikes structure...\n')
                obj.SpikeUMSStruct = [];
            end
        end
        
        function spksStruct = get.SpikeUMSStruct(obj)
            spksStruct = obj.SpikeUMSStruct;
        end
        
        %% GET SpikeTimes
        
        function spkTms = get.SpikeTimes(obj)
            spkTms = [];
            if ~isempty(obj.SpikeUMSStruct) &&...
                    isfield(obj.SpikeUMSStruct,'assigns')
                spikes = obj.SpikeUMSStruct;
            else
                disp(['The UMS2k pipeline hasn''t been ran or ',...
                    'the spike structure hasn''t been saved'])
                disp('No spike times extracted, therefore none returned')
                return;
            end
            % Return the spike times only if 1.- the cluster number exist;
            % 2.- the cluste is not garbage; or 3.- the cluster is a good
            % unit.
            lbls = obj.SpikeUMSStruct.labels;
            % goodSpks = lbls(:,2) == 2 | lbls(:,2) ~= 4;
            goodID = UMSDataLoader.getClustersID(obj.SpikeUMSStruct,'clType',...
                'good');
            if numel(goodID) == 1
                spkTms = spikes.spiketimes(spikes.assigns == goodID);
            elseif numel(goodID) > 1
                Ncls = 1;
                spkTms = cell(1,numel(goodID));
                for ccl = 1:numel(goodID)
                    spkTms(Ncls) =...
                        {spikes.spiketimes(spikes.assigns == goodID(ccl))};
                    Ncls = Ncls + 1;
                end
            else
                warning('All clusters were labeled as ''garbage''')
            end
            obj.SpikeTimes = spkTms;
        end
        
        function spkSubs = getSpikeSubscripts(obj)
            spkTms = obj.SpikeTimes;
            spkSubs = round(spkTms.*obj.SamplingFrequency);
        end
        
        %% SAVE and LOAD spiketimes
        function saveSpikeTimes(obj,currdir)
            if exist('currdir','var')
                [inDir,baseName,~] = fileparts(currdir);
            else
                [inDir, baseName, ~] = fileparts(obj.FileName);
            end
            sortFileName = [fullfile(inDir,baseName),'sorted.mat'];
            if ~isempty(obj.SpikeTimes)
                if ~exist(sortFileName,'file')
                    spikes = obj.SpikeUMSStruct;
                    spkTms = obj.SpikeTimes;
                    SPKS1 = struct('spikes',spikes,...
                        'spkTms',spkTms); %#ok<*NASGU>
                    try
                        save(sortFileName,'SPKS1')
                    catch 
                        disp('There was a problem saving the file')
                    end
                else
                    disp('The sorted file exists. We should discuss what to do in the lab meeting')
                    % The idea is to append a structure that contains the
                    % parameters used for the sorting and the extracted
                    % spike times. 
                end
            else
                return;
            end
        end
        
        function loadSpikeTimes(obj, fileName)
            if nargin == 2
                if exist(fileName,'file')
                    load(fileName,'SPKS1')
                    if isempty(obj.Data) && exist('SPKS1','var')
                        % It is recommended that the object is properly
                        % initialized before calling the importing method.
                        [inDir, baseName, fileExt] = fileparts(fileName);
                        dataBaseName = strsplit(baseName,'sorted');
                        dataFile = fullfile(inDir,[dataBaseName{1},fileExt]);
                        obj.changeDataFromFile(dataFile)
                    end
                else
                    fprintf('Wrong file name!\n')
                end
            elseif nargin == 1 && ~isempty(obj.FileName) &&...
                    ~isempty(obj.Data)
                [inDir,baseName,~] = fileparts(obj.FileName);
                sortedFile = fullfile(inDir,[baseName,'sorted.mat']);
                if exist(sortedFile,'file')
                    load(sortedFile,'SPKS1')
                else
                    fprintf('Bad news: the sorted file doesn''t exist\n')
                    return
                end
            else
                fprintf('This function accepts zero or one input argument.\n')
                return
            end
            fprintf('Loading UMS spikes structure and spike times...\n')
            obj.SpikeUMSStruct = SPKS1.spikes;
            obj.SpikeTimes = SPKS1.spkTms/obj.SamplingFrequency;
            fprintf('Done!\n')
        end
        
        %% Poly channel order modification
        function changeChannelOrder(obj,newChanOrd)
            obj.PolyChanOrder = newChanOrd;
            obj.ChannelsID = newChanOrd;
        end
        
        function set.PolyChanOrder(obj,newChanOrd)
            % The channel order will change if the Data property is empty,
            % if the user agrees to change the order even if the number of
            % channels
            if ~isempty(obj.Data) %#ok<*MCSUP>
                if obj.Nch == length(newChanOrd)
                    obj.PolyChanOrder = newChanOrd;
                else
                    fprintf(['WARNING! The data has %d channels and the',...
                        ' new order has %d numbers.'], obj.Nch,...
                        length(newChanOrd))
                    yn = input('Would you like to continue?','s');
                    if yn == 'y' || yn == 'Y' || strcmpi(yn,'yes')
                        obj.PolyChanOrder = newChanOrd;
                        disp(obj)
                    end
                end
            else
                obj.PolyChanOrder = newChanOrd;
            end
        end
        %% Get dependent variables
        % Number of samples per channel
        function samplesNumber = get.Ns(obj)
            samplesNumber = 0;
            if ~isempty(obj.Data)
                samplesNumber = size(obj.Data{1},1);
            end
        end
        
        % Number of channels in the data set
        function numberOfChannels = get.Nch(obj)
            numberOfChannels = 0;
            if ~isempty(obj.Data{1})
                numberOfChannels = size(obj.Data{1},2);
            end
        end
        
        function timeOut = get.Time(obj)
            if ~isempty(obj.Data)
                timeOut = (0:obj.Ns-1)/obj.SamplingFrequency;
            else
                timeOut = 0;
            end
        end
        
        %% Display of the object. Called everytime there is no semicolon.
        function disp(obj)
            fprintf('---UMSDataLoader object---\n')
            fprintf('File name: %s\n',obj.FileName)
            fprintf('Samples: %d\n',obj.Ns)
            fprintf('Sampling Frequency: %.3f kHz\n',obj.SamplingFrequency/1e3)
            fprintf('Detection method: %s\n',obj.DetectMethod)
            fprintf('Threshold: %f\n',obj.Thresh)
            if ~isempty(obj.Data)
                fprintf('Number of channels: %d\n',obj.Nch)
            end
            fprintf('Channel order: \n')
            for cch = 1:length(obj.PolyChanOrder)
                fprintf('Channel %d has ID %d from the file\n',cch,obj.PolyChanOrder(cch))
            end
            fprintf('Spike Structure: ')
            if ~isempty(obj.SpikeUMSStruct)
                fprintf('1\n')
            else
                fprintf('0\n')
            end
        end
        
        % Plot the data in the object on different levels
        function h = plot(obj,varargin)
            if ~isempty(obj.Data{1})
                h = figure('Name','UltraMegaSort2000 data','Color',[1,1,1]);
                means = mean(obj.Data{1},1);
                stds = std(obj.Data{1},[],1);
                if length(obj.PolyChanOrder) == obj.Nch
                    tx = obj.Time;
                    lbls = cell(1,obj.Nch);
                    lvl = cumsum(stds*30);
                    for cch = 1:obj.Nch
                        lbls(cch) = {sprintf([num2str(cch),' (ID ',...
                            num2str(obj.PolyChanOrder(cch)),...
                            ' ', char(obj.ChannelsID(cch)),')'])};
                        tempChan = obj.Data{1}(:,cch);
                        tempChan = tempChan -...
                            means(cch) + lvl(cch);
                        plot(tx,tempChan,varargin{:})
                        if cch == 1
                            hold on
                        end
                        NTcl = obj.SpikeUMSStruct.labels(:,1);
                        cmap = jet(numel(NTcl));
                        clID = UMSDataLoader.getClustersID(...
                            obj.SpikeUMSStruct,'clType','good');
                        [~,clSub] = intersect(NTcl,clID);
                        if isnumeric(obj.SpikeTimes)
                            spTms = obj.SpikeTimes;
                            plot(spTms, tempChan(...
                                round(spTms*obj.SamplingFrequency)),...
                                'LineStyle','none','Marker','^',...
                                'DisplayName',sprintf('Unique cluster: %d',...
                                clID),...
                                'Color',cmap(clSub,:))
                        elseif isa(obj.SpikeTimes,'cell')
                            for ccl = 1:numel(obj.SpikeTimes)
                                cSpTms = obj.SpikeTimes{ccl};
                                plot(cSpTms,tempChan(...
                                    round(cSpTms*obj.SamplingFrequency)),...
                                    'LineStyle','none','Marker','^',...
                                    'DisplayName',sprintf('Cluster %d',...
                                    clID(ccl)),...
                                    'Color',cmap(clSub(ccl),:))
                            end
                        end
                        legend('show','Location','best')
                    end
                    hold off;box off
                    set(gca,'YTick',lvl,'YTickLabel',lbls,...
                        'TickLabelInterpreter','none')
                    xlabel('Time [s]');title('Loaded data')
                else
                    fprintf(...
                        ['The number of channels in the data is ',...
                        'different from the order of channels!\nMaybe',...
                        ' change the Chanel Order property using the',...
                        ' changeChanOrder method.\nThe available ',...
                        'channels are the following:\n'])
                    for dispCount = 1:obj.Nch
                        if iscell(obj.ChannelsID)
                            fprintf('%s\n',obj.ChannelsID{dispCount})
                        else
                            fprintf('%d\n',obj.ChannelsID(dispCount))
                        end
                    end
                end
            end
        end
        
        % Plot the 'good-unit' waveforms
        function p = plotWaveforms(obj,clType)
            if isempty(obj.SpikeUMSStruct) || ~isfield(obj.SpikeUMSStruct,'labels')
                fprintf(1,'The UMS structure is empty!\n')
                fprintf(1,'Seems like the UMS pipeline has not been run\n')
                fprintf(1,'Type:\n')
                fprintf(1,'obj.UMS2kPipeline\n')
                fprintf(1,'to run the UMS2000 pipeline\n')
                return;
            end
            if exist('clType','var')
                clType = lower(clType);
            else
                clType = 'all';
            end
            try
                clID = UMSDataLoader.getClustersID(...
                    obj.SpikeUMSStruct,'clType',clType);
            catch ME
                fprintf(1,'It is likely that either the UMS2000 pipeline')
                fprintf(1,' or that the cluster type was not correct\n')
                fprintf(1,'Input given: %s\t',clType)
                fprintf(1,'Valid inputs:\n''good'',\n''mua'',\n')
                fprintf(1,'''unsorted'',\n''garbage''\n')
                disp(ME)
                return
            end
            Ncl = numel(clID);
            Nr = floor(sqrt(Ncl));
            Nc = ceil(Ncl/Nr);
            fs = obj.SamplingFrequency;
            Nsamp = size(obj.SpikeUMSStruct.waveforms,2);
            tx = (0:Nsamp-1)/fs;
            waveformFig = figure('Name',sprintf('%s clusters',clType),...
                'Visible','off','Color',[1,1,1]);
            ax = gobjects(1,Ncl);
            p = ax;
            
            NTcl = obj.SpikeUMSStruct.labels(:,1);
            cmap = jet(numel(NTcl));
            [~,clSub] = intersect(NTcl,clID);

            for cwf = 1:Ncl
                ax(cwf) = subplot(Nr,Nc,cwf,'Parent',waveformFig);
                clIdx = obj.SpikeUMSStruct.assigns == clID(cwf);
                cclWaveforms = obj.SpikeUMSStruct.waveforms(clIdx,:);
                plot(ax(cwf),tx,cclWaveforms,'LineWidth',0.2,...
                    'Color',[0.5,0.5,0.51]);
                ax(cwf).NextPlot = 'add';
                p(cwf) = plot(ax(cwf),tx,mean(cclWaveforms,1),...
                    'LineWidth',3,'Color', cmap(clSub(cwf),:),...
                    'DisplayName',sprintf('Cluster %d',clID(cwf)));
                title(ax(cwf),sprintf('Cluster %d',clID(cwf)));
                xlabel(ax(cwf),'Time [s]')
                ylabel(ax(cwf),'Amplitude [V]')
            end
            waveformFig.Visible = 'on';
        end
    end % methods
    %% Private methods
    methods (Access = 'private')
        % This function is a helper function to extract the spikes
        % structure out of the splitmerge_tool GUI from UMS.
        function getSpikeStructureCls(obj,varargin)
            try
                h = varargin{end};
                figdata = get(h,'UserData');
                obj.SpikeUMSStruct = figdata.spikes;
            catch
                disp('There was a problem reading the spike structure')
                disp('Please start an issue process in GitHub...')
                disp('Sorry :(')
            end
            delete(h)
        end
    end
    %% Static and private methods.
    methods (Static, Access = 'private')
        function [data, fs, chanReadOut, fileNameOut, chanTitle] =...
                LoadFromFile(fileName,channelOrder)
            % Private static method which reads the channels from the .mat
            % file created from the .smr or .smrx files with some metadata.
            try
                chanVars = load(fileName,'chan*');
                headVars = load(fileName,'head*');
            catch
                % If the file is not found, the given name is compared
                % against all the files in the given directory and in case
                % that the name is 'close' enough, it is presented to the
                % user.
                disp('No recognizable variables nor file!')
                [inDir,fname,ftype] = fileparts(fileName);
                filesInDir = dir(inDir);
                simFile = zeros(1,numel(filesInDir)-2);
                for cf = 3:numel(filesInDir)
                    simFile(cf-2) = sum(diag(...
                        distmatrix(...
                        (filesInDir(cf).name)',...
                        ([fname,ftype])')));
                end
                [~,shortDist] = min(simFile);
                disp('Only .mat files will be recognized!')
                errCorr = input(['Did you mean: ',filesInDir(shortDist+2).name,...
                    '? (y/n)'], 's');
                if errCorr == 'y' || errCorr == 'Y'
                    fileNameOut = fullfile(inDir,filesInDir(shortDist+2).name);
                    [data, fs, chanReadOut, fileNameOut] = UMSDataLoader.LoadFromFile(...
                        fileNameOut, channelOrder);
                else
                    disp('Loading aborted!')
                    data = [];
                    fs = 0;
                    chanReadOut = 1;
                end
                return
            end % Try to construct or compile the channels together.
            Ns = max(structfun(@numel,chanVars));
            szs = structfun(@length,chanVars);
            mu = mean(szs);
            sig = std(szs);
            if sig/max(szs) > 0.05
                signalFlag = szs > mu;
            else
                signalFlag = true(1,numel(szs));
            end
            chanNames = fieldnames(chanVars);
            headNames = fieldnames(headVars);
            chanIDs = arrayfun(@str2double,...
                cellfun(@(x) x(5:end),chanNames,'UniformOutput',false));
            if sum(signalFlag) ~= numel(signalFlag)
                fprintf('Not all channels in file are continuous signals!\n')
                fprintf('Channel(s) ')
                for cch = find(~signalFlag)
                    fprintf('%d ',chanIDs(cch))
                end
                fprintf('is (are) not signals.\n')
            end
            fprintf('Deleting...\n')
            chanIDs = chanIDs(signalFlag);
            chanTitle = cell(1,numel(chanIDs));
            if ~exist('channelOrder','var') || isempty(channelOrder)
                channelOrder = chanIDs;
            end
            dataMatrix = zeros(Ns,numel(channelOrder));
            chanReadOut = zeros(1,numel(chanIDs));
            chanKeepIdx = true(1,numel(channelOrder));
            Nsamples = zeros(1,numel(chanIDs));
            for cch = 1:numel(channelOrder)
                inFl = chanIDs == channelOrder(cch);
                % Something important to consider is that the importing
                % could be taking repeated channels if and only if they
                % exist in the file.
                if sum(inFl)
                    auxChan = chanVars.(chanNames{inFl});
                    Nsamples(cch) = length(auxChan);
                    auxHead = headVars.(headNames{inFl});
                    chanTitle(cch) = {auxHead.title};
                    if isrow(auxChan)
                        auxChan = auxChan';
                    end
                    if cch == 1
                        try
                            fs = auxHead.SamplingFrequency;
                        catch
                            disp(['No information about the ',...
                                'sampling frequency in this file.'])
                            fs = 2e4;
                            fprintf('Setting sampling frequency: %0.1f kHz\n',...
                                fs/1e3)
                        end
                    end
                    chanReadOut(cch) = channelOrder(cch);
                    dataMatrix(1:Nsamples(cch),cch) = auxChan;
                else
                    chanKeepIdx(cch) = false;
                    fprintf('Channel %d is not part of this recording!\n',...
                        channelOrder(cch))
                end
            end
            Nsamples(Nsamples==0) = [];
            Ns = min(Nsamples);
            chanTitle(chanReadOut==0) = [];
            chanReadOut(chanReadOut==0) = [];
            fileNameOut = fileName;
            data{1} = dataMatrix(1:Ns,chanKeepIdx);
        end
    end % Static methods
    %% Static and public methods
    methods (Static, Access = 'public')
        function outID = getClustersID(spkUMSstr, varargin)
            % Parsing inputs
            outID = [];
            p = inputParser;
            p.CaseSensitive = true;
            p.KeepUnmatched = true;
            defaultTypes = 'all';
            validTypes = {'good','mua','unsorted','garbage','all'};
            checkTypes = @(x) any(validatestring(x,validTypes));
            
            validStruct = @(x) isstruct(x) && isfield(x,'labels');
            
            addRequired(p,'spkUMSstr',validStruct)
            addOptional(p,'clType',defaultTypes,checkTypes)
            try
                parse(p,spkUMSstr,varargin{:});
            catch EM
                fprintf(1,'There was something wrong with the inputs...\n')
                disp(EM)
                return
            end
            spkUMSstr = p.Results.spkUMSstr;
            clType = p.Results.clType;
            p.delete;
            clearvars p
            clID = spkUMSstr.labels(:,1);
            switch clType
                case 'all'
                    outID = clID(spkUMSstr.labels(:,2) ~= 4);
                case 'good'
                    outID = clID(spkUMSstr.labels(:,2) == 2);
                case 'mua'
                    outID = clID(spkUMSstr.labels(:,2) == 3);
                case 'unsorted'
                    outID = clID(spkUMSstr.labels(:,2) == 1);
                case 'garbage'
                    outID = clID(spkUMSstr.labels(:,2) == 4);
            end
        end
    end
end % classdef


