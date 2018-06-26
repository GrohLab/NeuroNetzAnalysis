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
        Shadow (1,1) double= 0.9;            % ms, enforced dead region after each spike
        RefractoryPeriod (1,1) double= 2.5;  % ms, refractory period (for calculation refractory period violations)
        % Alignment parameters
        WindowSize (1,1) double = 1.5;       % ms, width of a spike
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
        
        function invertSignals(obj)
            obj.Data{1} = -obj.Data{1};
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
                obj.SpikeTimes = [];
                obj.SpikeUMSStruct = [];
            end
        end
        
        %% SET AND GET SpikeUMSStruct
        function set.SpikeUMSStruct(obj,structIn)
            if isstruct(structIn) && isfield(structIn,'params')
                obj.SpikeUMSStruct = structIn;
            else
                fprintf('Deleting the UMS spikes structure...\n')
                obj.deleteSpikeUMSStruct;
            end
        end
        
        function spksStruct = get.SpikeUMSStruct(obj)
            spksStruct = obj.SpikeUMSStruct;
        end
        
        function deleteSpikeUMSStruct(obj)
            obj.SpikeUMSStruct = [];
        end
        
        %% GET, SET and SAVE SpikeTimes
        function spksTime = get.SpikeTimes(obj)
            %#ok<*MCSUP>
            if ~isempty(obj.SpikeTimes)
                spksTime = obj.SpikeTimes;
            else
                spksTime = [];
                % disp('No spike times extracted yet');
            end
        end % get SpikeTimes
        
        function getSpikeTimes(obj)
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
            goodSpks = lbls(:,2) == 2 | lbls(:,2) ~= 4;
            
            obj.SpikeTimes =...
                spikes.spiketimes(spikes.assigns == lbls(goodSpks,1));
        end
        
        function saveSpikeTimes(obj)
            [inDir, baseName, ~] = fileparts(obj.FileName);
            sortFileName = [fullfile(inDir,baseName),'sorted.mat'];
            if ~isempty(obj.SpikeTimes)
                if ~exist(sortFileName,'file')
                    spikes = obj.SpikeUMSStruct;
                    spkTms = obj.SpikeTimes;
                    SPKS1 = struct('spikes',spikes,...
                        'spkTms',round(obj.SamplingFrequency*spkTms));
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
            if ~isempty(obj.Data{1})
                timeOut = seconds(0:1/obj.SamplingFrequency:...
                    (obj.Ns-1)/obj.SamplingFrequency);
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
                        lbls{cch} = [num2str(cch),' (ID ',...
                            num2str(obj.PolyChanOrder(cch)),')'];
                        tempChan = obj.Data{1}(:,cch);
                        tempChan = tempChan -...
                            means(cch) + lvl(cch);
                        plot(tx,tempChan,varargin{:})
                        if cch == 1
                            hold on
                        end
                        if ~isempty(obj.SpikeTimes)
                            spIx = round(obj.SamplingFrequency*obj.SpikeTimes);
                            plot(tx(spIx), tempChan(spIx),...
                                'LineStyle','none','Marker','.')
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
                h = varargin{end};
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
            Ns = max(structfun(@numel,chanVars));
            chanNames = fieldnames(chanVars);
            headNames = fieldnames(headVars);
            chanIDs = arrayfun(@str2double,...
                cellfun(@(x) x(5:end),chanNames,'UniformOutput',false));
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
                    dataMatrix(:,cch) = auxChan;
                else
                    chanKeepIdx(cch) = false;
                    fprintf('Channel %d is not part of this recording!\n',...
                        channelOrder(cch))
                end
            end
            Nsamples(Nsamples==0) = [];
            Ns = min(Nsamples);
            chanReadOut(chanReadOut==0) = [];
            fileNameOut = fileName;
            data{1} = dataMatrix(1:Ns,chanKeepIdx);
        end
    end % Static methods
end % classdef


