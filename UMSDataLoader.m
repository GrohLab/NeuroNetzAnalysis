classdef UMSDataLoader
    % DATALOADER 
    properties
        Data = [];
    end % properties
    properties (SetAccess = 'private',GetAccess = 'private')
        % Default value is for the Poly design unaccesible to the user:
        PolyChanOrder = [8, 9, 7, 10, 4, 13, 5, 12, 2, 15, 1, 16, 6, 11, 3, 14];
    end
    methods
        function obj = UMSDataLoader(varargin)
            if nargin == 1
                if isnumeric(varargin{1})
                    if iscolumn(varargin{1})
                        data = double(varargin{1})';
                    else
                        data = double(varargin{1});
                    end
                elseif ischar(varargin{1})
                    disp('Building from file...')
                    auxObj = UMSDataLoader();
                    data = UMSDataLoader.LoadFromFile(...
                        varargin{1},auxObj.PolyChanOrder);
                else
                    error('Sorry! I wouldn''t know how to deal with this input')
                end
            elseif nargin >= 2
                disp('Preparing this code')
                data = [];
            else
                data = [];
            end
            obj.Data = data;
        end
        function spikes = UMS2000(obj,samplingFreq,thresh)
            if nargin == 2
                spikes = ss_default_params(samplingFreq);
            elseif nargin == 3
                spikes = ss_default_params(samplingFreq,thresh);
            end
            spikes = ss_detect(obj.Data,spikes);
            spikes = ss_align(spikes);
            spikes = ss_kmeans(spikes);
            spikes = ss_energy(spikes);
            spikes = ss_aggregate(spikes);
            splitmerge_tool(spikes)
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
                end
            end
            
        end
    end
end % classdef

% data{1} = fR';
% spikes = ss_default_params(fs);
% spikes = ss_detect(data,spikes);
% spikes = ss_align(spikes);
% spikes = ss_kmeans(spikes);
% spikes = ss_energy(spikes);
% spikes = ss_aggregate(spikes);
