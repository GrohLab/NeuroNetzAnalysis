classdef SMRXConverter < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = 'protected')
        SmrxFiles
        BinFiles
    end
    
    methods
        function obj = SMRXConverter(DataDir)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            binFiles = dir(fullfile(DataDir, '*.bin'));
            smrxFiles = dir(fullfile(DataDir, '*.smrx'));
            if isempty(smrxFiles)
                fprintf(1,'Seems like this folder does not correspond')
                fprintf(1,' to an experiment folder.\n')
                fprintf(1,'Please be sure that all experiment relevant')
                fprintf(1,' files are in the chosen folder.\n')
                obj.delete;
                return
            end
            if numel(smrxFiles) > 1
                fprintf(1,'More than one .smrx files are in this folder.\n')
                str = input(['Would you like to open a GUI to ',...
                    'select one (or many to concatenate)? (y/n)'],'s');
                if ~isItYes(str)
                    fprintf(1,'No .bin file created!\n')
                    return
                end
                smrxNames = arrayfun(@(x) x.name, smrxFiles,...
                    'UniformOutput', 0);
                [chSmrx, iOk] = listdlg('ListString', smrxNames,...
                    'SelectionMode', 'multiple');
                if ~iOk
                    fprintf(1,'Cancelled.\n');
                    return
                end
                Nf = length(chSmrx);
                if Nf > 1
                    fileOrder = (1:Nf)';
                    defInput = num2cell(num2str(fileOrder));
                    answr = inputdlg(smrxNames,'File order',[1, 60],defInput);
                    nFileOrder = str2double(answr);
                    nSmrxFiles = smrxFiles;
                    if ~isempty(answr) && sum(abs(fileOrder - nFileOrder)) ~= 0
                        fprintf(1,'Changing file order...\n')
                        nSmrxFiles(nFileOrder) = smrxFiles(fileOrder);
                        smrxFiles = nSmrxFiles;
                    else
                        fprintf('File order not altered\n')
                    end
                end
                smrxPaths = arrayfun(@(x) fullfile(x.folder, x.name),...
                    smrxFiles, 'UniformOutput', 0);
                
            end
            
            if isempty(binFiles)
            end
        end
        
        function outputArg = createBinFile(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
        end
    end
end

