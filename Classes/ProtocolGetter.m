classdef ProtocolGetter
    %PROTOCOLGETTER tries to recognize and extract the experimental
    %protocol
    
    properties
        dataDir char
        ismerged logical
        fileOrder cell
    end
    
    methods
        function obj = ProtocolGetter(directory)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            if ~exist(directory,'dir')
                fprintf('Invalid folder! No object created!\n')
                return
            end
            obj.dataDir = directory;
            binFiles = dir([directory, '*.bin']);
            smrxFiles = dir([directory, '*.smrx']);
            Nb = numel(binFiles); Ns = numel(smrxFiles);
            [~,binBaseNames] = arrayfun(@(x) fileparts(x.name), binFiles,...
                'UniformOutput', 0);
            [~,smrxBaseNames] = arrayfun(@(x) fileparts(x.name), smrxFiles,...
                'UniformOutput', 0);
            % Checking if the number of binnary files equals the number of
            % CED recorded files and if their names match. sitFlags has 4
            % likely cases encoded from
            sitFlags = [Nb == Ns,...
                any(ismember(binBaseNames, smrxBaseNames))];
            selNum = bin2num(sitFlags,0);
            switch selNum
                case 0 % Merged binnary file using all SMRX files; ideal
                    merFiles = dir([directory,'*_fileOrder.mat']);
                    if ~isempty(merFiles)
                        
                    end
                case 1 
                case 2
                case 3
            end
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

