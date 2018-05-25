classdef UMSDataLoader < handle
    % DATALOADER 
    properties (SetAccess = 'private')
        Data = [];
        SpikeTimes
    end
    properties
        SpikeUMSStruct
    end
    properties
        SamplingFrequency = 2e4;
        % Detecting parameters
        DetectMethod = 'auto';
        Thresh = 3.9;
        Shadow = 0.75;            % ms, enforced dead region after each spike
        RefractoryPeriod = 2.5;  % ms, refractory period (for calculation refractory period violations)
        % Alignment parameters
        WindowSize  = 1.5;       % ms, width of a spike
        CrossTime = 0.6;         % ms, alignment point for peak of waveform
        MaxJitter = 0.6;         % ms, width of window used to detect peak after threshold crossing
        % sorting parameters
        AggCutoff = .05;         %  higher = less aggregation, lower = more aggregation
        KmeansClusterSize = 500; %  target size for miniclusters
    end % properties
    properties (Dependent)
        Ns 
    end
    properties (SetAccess = 'private',GetAccess = 'private')
        % Default value is for the Poly design unaccesible to the user:
        PolyChanOrder = [8, 9, 7, 10, 4, 13, 5, 12, 2, 15, 1, 16, 6, 11, 3, 14];
    end
    methods
        function obj = UMSDataLoader(varargin)
            % For UMS2k to process the spike traces, it needs a cell array
            % of Nt (number of trials) which contains a Ns x Nch (number of
            % samples x number of channels) matrix. This object arranges
            % the given data into such format from a file or from the
            % workspace.
            if nargin == 1
                if isnumeric(varargin{1})
                    if isrow(varargin{1})
                        data{1} = double(varargin{1})';
                    else
                        data{1} = double(varargin{1});
                    end
                elseif ischar(varargin{1})
                    disp('Building from file...')
                    data = UMSDataLoader.LoadFromFile(...
                        varargin{1},obj.PolyChanOrder);
                elseif iscell(varargin{1})
                   data = varargin{1};
                else
                    error('Sorry! I wouldn''t know how to deal with this input')
                end
            elseif nargin >= 2
                disp('Joining the channels into a matrix...')
                disp('(min_chan_length x #_channels)')
                try 
                    dataMat = double(cell2mat(varargin));
                    [Ns, Nch] = size(dataMat);
                    if Ns < Nch
                        dataMat = dataMat';
                    end
                    data{1} = dataMat;
                catch
                    szs = cell2mat(cellfun(@size,varargin,...
                        'UniformOutput',false)');
                    Ns = min(szs(szs ~= 1));
                    dataMatrix = zeros(Ns,nargin);
                    for carg = 1:nargin
                        temp = varargin{carg};
                        if iscolumn(temp)
                            temp = temp';
                        end
                        dataMatrix(:,carg) =  temp;
                    end
                    data{1} = double(dataMatrix);
                end
            else
                data = [];
            end
            obj.Data = data;
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
            spikesLocal = ss_detect(obj.Data,spikesLocal);
            spikesLocal = ss_align(spikesLocal);
            spikesLocal = ss_kmeans(spikesLocal);
            spikesLocal = ss_energy(spikesLocal);
            obj.SpikeUMSStruct = ss_aggregate(spikesLocal);
            splitmerge_tool(obj.SpikeUMSStruct)
            
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
        
        function obj = getSpikeTimes(obj,varargin)
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
        
        function samplesNumber = get.Ns(obj)
            samplesNumber = numel(obj.Data{1});
        end
        function disp(obj)
            fprintf('UMSDataLoader object:\nSamples: %d\n',obj.Ns)
            fprintf('Sampling Frequency: %.3f kHz\n',obj.SamplingFrequency/1e3)
            fprintf('Detection method: %s\n',obj.DetectMethod)
            fprintf('Threshold: %f\n',obj.Thresh)
        end
    end % methods
    
    methods (Static)
        function data = LoadFromFile(fileName,channelOrder)
            try
                disp('Loading precompiled data array...')
                load(fileName,'Data')
                data = Data;
            catch
                disp('Not found.')
                disp('Compiling channels...')
                try
                    chanVars = load(fileName,...
                        'chan*');
                    Ns = min(structfun(@numel,chanVars));
                    chanNames = fieldnames(chanVars);
                    dataMatrix = zeros(Ns,numel(channelOrder));
                    for cch = 1:numel(channelOrder)
                        auxChan = chanVars.(chanNames{channelOrder(cch)})(1:Ns);
                        if isrow(auxChan)
                            auxChan = auxChan';
                        end
                        dataMatrix(:,cch) = auxChan;
                    end
                    data{1} = dataMatrix;
                catch
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
                        auxObj = UMSDataLoader();
                        data = UMSDataLoader.LoadFromFile(...
                            fullfile(inDir,filesInDir(shortDist+2).name),...
                            auxObj.PolyChanOrder);
                    else
                        disp('Loading aborted!')
                        data = [];
                    end
                end % Try to construct or compile the channels together.
            end % Try to load the data cell from file
            
        end
    end % Static methods
end % classdef
