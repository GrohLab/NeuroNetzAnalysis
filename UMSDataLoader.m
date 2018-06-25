classdef UMSDataLoader < handle
    % DATALOADER
    properties (SetAccess = 'private')
        Data(1,:) cell = {};
        SpikeTimes
    end
    properties
        SpikeUMSStruct;
    end
    properties
        SamplingFrequency (1,1) double = 2e4;
        % Detecting parameters
        DetectMethod (1,1) char = 'auto';
        Thresh (1,1) single= 3.9;
        Shadow (1,1) single= 0.85;            % ms, enforced dead region after each spike
        RefractoryPeriod (1,1) single = 2.5;  % ms, refractory period (for calculation refractory period violations)
        % Alignment parameters
        WindowSize (1,1) single = 1.5;       % ms, width of a spike
        CrossTime (1,1) single = 0.6;         % ms, alignment point for peak of waveform
        MaxJitter (1,1) single = 0.6;         % ms, width of window used to detect peak after threshold crossing
        % sorting parameters
        AggCutoff (1,1) single = .05;         %  higher = less aggregation, lower = more aggregation
        KmeansClusterSize (1,1) int8 = 500; %  target size for miniclusters
    end % properties
    properties (Dependent)
        Ns                      % Number of samples
        Nch                     % Number of channels
    end
    properties (SetAccess = 'private',GetAccess = 'private')
        % Default value is for the Poly design unaccesible to the user:
        PolyChanOrder (1,:) int16 =...
            [8, 9, 7, 10, 4, 13, 5, 12, 2, 15, 1, 16, 6, 11, 3, 14];
        FileName char;
    end
    methods
        function obj = UMSDataLoader(filename,chanOrder)
            % For UMS2k to process the spike traces, it needs a cell array
            % of Nt (number of trials) which contains a Ns x Nch (number of
            % samples x number of channels) matrix. This object arranges
            % the given data into such format from a file or from the
            % workspace.
            if nargin == 1
                [data, fs, chanReadOut, fileNameOut] = UMSDataLoader.LoadFromFile(...
                    filename,[]);
            elseif nargin == 2
                [data, fs, chanReadOut, fileNameOut] = UMSDataLoader.LoadFromFile(...
                    filename,chanOrder);
            end
            if ~isempty(data)
                obj.FileName = fileNameOut;
                obj.Data = data;
                obj.changeChannelOrder(chanReadOut);
                obj.SamplingFrequency = fs;
            end
            
            disp('Constructed!')
        end % Constructor
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
            set(h,'CloseRequestFcn',{@obj.getSpikeStructureCls,h,strIdx})
            
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
                [data, fs, chanReadOut, fileNameOut] =...
                    obj.LoadFromFile(fileName,[]);
            elseif nargin == 3
                [data, fs, chanReadOut, fileNameOut] =...
                    obj.LoadFromFile(fileName,chanOrder); 
            end
            if ~isempty(data)
                obj.Data = data;
                obj.SamplingFrequency = fs;
                obj.changeChannelOrder(chanReadOut);
                obj.FileName = fileNameOut;
            end
        end
        
        
        
        %% SET AND GET SpikeUMSStruct
        function set.SpikeUMSStruct(obj,structIn)
            if isstruct(structIn) && isfield(structIn,'params')
                if ~isempty(obj.SpikeUMSStruct)
                    yon = input(['The structure is not empty, ',...
                        'do you wish to overwrite? (y/n)'],'s');
                    if yon == 'y' || yon == 'Y'
                        obj.SpikeUMSStruct = structIn;
                    end
                else
                    obj.SpikeUMSStruct = structIn;
                end
            else
                fprintf('The given input is not a valid structure.')
            end
        end
        
        function spksStruct = get.SpikeUMSStruct(obj)
            spksStruct = obj.SpikeUMSStruct;
        end
        %% GET & SET SpikeTimes
        function spksTime = get.SpikeTimes(obj)
            %#ok<*MCSUP>
            if ~isempty(obj.SpikeTimes)
                spksTime = obj.SpikeTimes;
            else
                spksTime = [];
                disp('No spike times extracted yet');
            end
        end % get SpikeTimes
        
        function getSpikeTimes(obj,varargin)
            spkClust = 0;
            for na = 1:nargin-1
                if isa(varargin{na},'struct')
                    spikes = varargin{na};
                elseif isa(varargin{na},'double')
                    spkClust = varargin{na};
                else
                    disp('Input ignored...')
                    disp(varargin{na})
                end
            end
            if exist('spikes','var') && isfield(spikes,'assigns')
                clst = unique(spikes.assigns);
            elseif ~isempty(obj.SpikeUMSStruct) &&...
                    isfield(obj.SpikeUMSStruct,'assigns')
                clst = unique(obj.SpikeUMSStruct.assigns);
                spikes = obj.SpikeUMSStruct;
            else
                disp(['The UMS2k pipeline hasn''t been ran or ',...
                    'the spike structure hasn''t been saved'])
                disp('No spike times extracted, therefore none returned')
                return;
            end
            if sum(clst == spkClust)
                obj.SpikeTimes =...
                    spikes.spiketimes(spikes.assigns == spkClust);
            else
                disp('The seleted cluster is non existing!')
                disp(['Please select a valid spike cluster. ',...
                    'HINT! You can see the cluster number in the splitmerge tool ;)'])
            end
        end
        %% Poly channel order modification
        function changeChannelOrder(obj,newChanOrd)
            obj.PolyChanOrder = newChanOrd;
        end
        
        function set.PolyChanOrder(obj,newChanOrd)
            % The channel order will change if the Data property is empty,
            % if the user agrees to change the order even if the number of
            % channels
            if ~isempty(obj.Data)
                if obj.Nch == length(newChanOrd)
                    obj.PolyChanOrder = newChanOrd;
                else
                    fprintf(['WARNING! The data has %d channels and the',...
                        ' new order has %d numbers.'], obj.Nch,...
                        length(newChanOrd))
                    yn = input('Would you like to continue?','s');
                    if yn == 'y' || yn == 'Y'
                        obj.PolyChanOrder = newChanOrd;
                        disp(obj)
                    end
                end
            else
                obj.PolyChanOrder = newChanOrd;
            end
        end
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
        
        % Display of the object. Called everytime there is no semicolon.
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
                fprintf('%d: %d\n',cch,obj.PolyChanOrder(cch))
            end
        end
        
        % Plot the data in the object on different levels
        function h = plot(obj,varargin)
            if ~isempty(obj.Data{1})
                h = figure('Name','UltraMegaSort2000 data','Color',[1,1,1]);
                means = mean(obj.Data{1},1);
                stds = std(obj.Data{1}(:,obj.PolyChanOrder),[],1);
                if length(obj.PolyChanOrder) >= obj.Nch
                    tx = seconds(0:1/obj.SamplingFrequency:...
                        (obj.Ns-1)/obj.SamplingFrequency);
                    lbls = cell(1,obj.Nch);
                    lvl = cumsum(stds(obj.PolyChanOrder)*30);
                    for cch = 1:obj.Nch
                        lbls{cch} = [num2str(cch),' (',...
                            num2str(obj.PolyChanOrder(cch)),')'];
                        tempChan = obj.Data{1}(:,obj.PolyChanOrder(cch));
                        tempChan = tempChan -...
                            means(obj.PolyChanOrder(cch)) + lvl(cch);
                        plot(tx,tempChan,varargin{:})
                        if cch == 1
                            hold on
                        end
                    end
                    hold off;box off
                    set(gca,'YTick',lvl,'YTickLabel',lbls)
                    xlabel('Time [s]');title('Loaded data')
                else
                    fprintf(...
                        ['The number of channels in the data is ',...
                        'different from the order of channels!\nMaybe',...
                        ' change the Chanel Order property using the',...
                        ' changeChanOrder method.'])
                end
            end
        end
    end % methods
    methods (Access = 'private')
        function getSpikeStructureCls(obj,varargin)
            try
                h = varargin{3};
                figdata = get(h,'UserData');
                obj.SpikeUMSStruct = figdata.spikes;
            catch
                disp('There was a problem reading the spike structure')
                disp('I am on debugging process... :( ')
            end
            delete(h)
        end
    end
    methods (Static)
        function [data, fs, chanReadOut, fileNameOut] =...
                LoadFromFile(fileName,channelOrder)
            try
                chanVars = load(fileName,...
                    'chan*');
                headVars = load(fileName,...
                    'head*');
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
            Ns = min(structfun(@numel,chanVars));
            chanNames = fieldnames(chanVars);
            headNames = fieldnames(headVars);
            chanIDs = arrayfun(@str2double,...
                cellfun(@(x) x(5:end),chanNames,'UniformOutput',false));
            chanTitle = cell(1,numel(chanIDs));
            if ~exist('channelOrder','var') || isempty(channelOrder)
                channelOrder = chanIDs;
            end
            load(fileName,['head',num2str(chanIDs(1))])
            dataMatrix = zeros(Ns,numel(channelOrder));
            chanReadOut = zeros(1,numel(chanIDs));
            chanKeepIdx = true(1,numel(channelOrder));
            for cch = 1:numel(channelOrder)
                inFl = chanIDs == channelOrder(cch);
                % Something important to consider is that the importing
                % could be taking repeated channels if and only if they
                % exist in the file.
                if sum(inFl)
                    auxChan = chanVars.(chanNames{chanIDs(inFl)})(1:Ns);
                    auxHead = headVars.(headNames{chanIDs(inFl)});
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
                    dataMatrix(:,cch) = auxChan;
                else
                    chanKeepIdx(cch) = false;
                    fprintf('Channel %d is not part of this recording!\n',...
                        channelOrder(cch))
                end
            end
            chanReadOut(chanReadOut==0) = [];
            fileNameOut = fileName;
            data{1} = dataMatrix(:,chanKeepIdx);
        end
    end % Static methods
end % classdef


