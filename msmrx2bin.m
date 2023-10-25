function iOk = msmrx2bin(dataDir, outBaseName, chanGroup)
iOk = -1;
impStr = 'Rhd';
istxt = @(x) ischar(x) | isstring(x);

if ~exist('chanGroup','var') || isempty(chanGroup)
    chanGroup = "";
end

smrxFiles = dir(fullfile(dataDir,'*.smrx'));
if isempty(smrxFiles)
    smrxFiles = dir(fullfile(dataDir,'*\*.smrx'));
    if isempty(smrxFiles)
        fprintf(1,'The given directory contains no .smrx files. Please ')
        fprintf(1,'try again with another folder which do contain .smrx files.\n')
        return
    end
end

cellSmrxFiles = {smrxFiles.name};
Nf = size(smrxFiles,1);

% Selecting files
[incFiles, iok] = listdlg('ListString',cellSmrxFiles(1,:),...
    'SelectionMode','multiple',...
    'PromptString','Select the files to join:',...
    'InitialValue',1:Nf);
if iok
    smrxFiles = smrxFiles(incFiles);
    cellSmrxFiles = {smrxFiles.name};
    Nf = size(smrxFiles,1);
else
    fprintf(1,'Cancelling...\n')
    return
end

% Merging order
fileOrder = (1:Nf)';
defInput = num2cell(num2str(fileOrder));
answr = inputdlg(cellSmrxFiles,'File order',[1, 60],defInput);
nFileOrder = str2double(answr);
nSmrxFiles = smrxFiles;
if ~isempty(answr) && sum(abs(fileOrder - nFileOrder)) ~= 0
    fprintf(1,'Changing file order...\n')
    nSmrxFiles(nFileOrder) = smrxFiles(fileOrder);
    smrxFiles = nSmrxFiles;
else
    fprintf('File order not altered\n')
end
clearvars nSmrxFiles nFileOrder
% Creating a .bin file
if ~exist('outBaseName','var') || isempty(outBaseName) ||...
        ~istxt(outBaseName)
    fprintf(1,'No outname given. Computing a name...\n')
    pathPieces = strsplit(dataDir,filesep);
    outBaseName = [pathPieces{end},'.bin'];
    fprintf(1,'File name: %s.bin\n',pathPieces{end})
end

outFullName = fullfile(dataDir,outBaseName);
if exist(outFullName,'file')
    ovwtAns = questdlg(...
        sprintf('Warning! File %s exists. Overwrite?',outFullName),...
        'Overwrite?','Yes','No','No');
    if strcmpi(ovwtAns,'no')
        fprintf('No file written.\n')
        return
    end
end
fID = fopen(outFullName,'w');
m = (2^32)/100;
fs = zeros(Nf,1);
for cf = 1:Nf
    cfName = fullfile(smrxFiles(cf).folder,smrxFiles(cf).name);
    FileInfo = SONXFileHeader(cfName);
    fhand = CEDS64Open(cfName);
    if fhand < 0
        fprintf(1,'The file might be opened in Spike2. Please close it and')        
        fprintf(1,' try again.\n')
        fclose(fID);
        return
    end
    totalTime = CEDS64TicksToSecs(fhand,FileInfo.maxFTime);
    
    if fhand > 0
        % Import the data.
        mxChans = CEDS64MaxChan(fhand);
        chanList = 1:mxChans;
%         try
%             heads = SONXChannelInfo(fhand,1,1);
%             fch = 2;
%             ch1 = 2;
%             chTypes(1) = 1;
%         catch
%             heads = SONXChannelInfo(fhand,2,2);
%             fch = 3;
%             ch1 = 3;
%             chTypes(2) = 1;
%         end
%         heads = repmat(heads,mxChans,1);
        chTypes = arrayfun(@(x) CEDS64ChanType(fhand, x), chanList');
        heads = arrayfun(@(x) SONXChannelInfo(fhand,x,x), chanList',...
            'UniformOutput', false);
        heads(chTypes~=1) = []; heads = cat(1, heads{:});
        heads = arrayfun(@(x) setfield(heads(x), "FileChannel", x),...
            (1:size(heads,1))');
%         for ch = ch1:mxChans
%             if chTypes(ch) == 1
%                 fch = fch + 1;
%                 try
%                     heads(ch) = SONXChannelInfo(fhand,ch,fch);
%                 catch
%                     auxHead = SONXChannelInfo(fhand,ch,fch);
%                     hdsFN = fieldnames(heads); auxFN = fieldnames(auxHead);
%                     [inH, inA] = ismember(hdsFN, auxFN);
% 
%                 end
%             end
%         end
        chanNames = string({heads.title}');
        desChans = contains(chanNames, chanGroup);
        heads = heads(desChans);
        chanList = chanList(chTypes == 1); chanList = chanList(desChans);
        chead = numel(heads);
        while chead >= 1
            if ~xor(isnan(str2double(heads(chead).title)),...
                    ~contains(heads(chead).title,impStr))
                heads(chead) = [];
                chanList(chead) = [];
            end
            chead = chead - 1;
        end
        multiplexerFactor = heads(1).ChanDiv;
        fs(cf) = 1 / (FileInfo.usPerTime * multiplexerFactor);
        FileInfo.SamplingFrequency = fs(cf);
        display(FileInfo)
        % Determining the necessary array size to occupy approximately the
        % 75% of the available memory given that the array is int16
        memStruct = memory;
        BuffSize = 3 * memStruct.MemAvailableAllArrays / 4;
        dataPointsExp = (BuffSize / (numel(chanList) * 2));
        if heads(1).npoints < dataPointsExp
            dataPointsExp = heads(1).npoints;
        end
        wwidth = double(dataPointsExp)/fs;
        % dataPointsExp = ceil(log10(fs)+2);
        % wwidth = 10^ceil(log10(fs)+2)/fs;
        is = 1./fs;
        cw = 0;
        if abs(totalTime-(heads(1).stop - heads(1).start))
            oldTotalTime = totalTime;
            totalTime = heads(1).stop - heads(1).start;
            fprintf(1,'The length of the signals in the file seem to ')
            fprintf(1,'differ (%.3f s \\delta (%.3f ms)).\nConsidering %.3f seconds\n',...
                oldTotalTime, 1e3*(oldTotalTime - totalTime), totalTime)
        end
        while cw < totalTime
            %         if exist('Npts','var')
            %             fID = fopen(outfilename,'a');
            %         end
            % dataBuff = zeros(numel(chanList),10^dataPointsExp,'int16');
            dataBuff = zeros(numel(chanList), dataPointsExp, 'int16');
            if cw <= totalTime - wwidth
                timeSegment = [cw, cw + wwidth];
            else
                timeSegment = [cw, totalTime];
                cw = totalTime * 2;
                dataBuff = zeros(numel(chanList),int32(diff(timeSegment)*fs(cf)),...
                    'int16');
            end
            shortFlag = false;
            fprintf(1,'Reading... ')
            for ch = 1:numel(chanList)
                [Npts, chanAux, ~] =...
                    SONXGetWaveformChannelSegment(fhand, chanList(ch), timeSegment,...
                    heads(ch)); %#ok<ASGLU>
                dat = int16(chanAux * m);
                try
                    dataBuff(ch,:) = dat;
                catch
                    dataBuff(ch,1:length(dat)) = dat;
                    shortFlag = true;
                end
            end
            
            fprintf(1,'done!\n')
            cw = cw + wwidth + is;
            fprintf(1,'Writting... ')
            if shortFlag
                dataBuff(:,length(dat)+1:dataPointsExp) = [];
            end
            % Removing the median of all channels.
            % TODO: 1.- save median per experiment
            % 2.- ask user for removing median.
            buffMedian = median(dataBuff);
            dataBuff = dataBuff - buffMedian;
            fwrite(fID,dataBuff,'int16');
            fprintf(1,' done!\n')
            %         ftell(fID)
            %         fclose(fID);
        end
    end
    CEDS64Close(fhand);
end
fclose(fID);
fs = mean(fs); 
[~, outBaseName] = fileparts(outFullName);
save(fullfile(dataDir,...
            string(outBaseName) + "_sampling_frequency.mat"),'fs')
fprintf(1, 'Successfully imported files!\n')
fprintf(1, 'The files merged are the following:\n')
smrxFileNames = cell(numel(smrxFiles),1);
ffoID = fopen(fullfile(dataDir,string(outBaseName) + "_fileOrder.txt"),'w');
for cf = 1:numel(smrxFiles)
    fprintf(1,'%s\n',smrxFiles(cf).name);
    smrxFileNames(cf) = {smrxFiles(cf).name};
    fprintf(ffoID, '%s\n', smrxFiles(cf).name);
end
fprintf(ffoID, '%s.bin', outBaseName);
fclose(ffoID);
iOk = CEDS64Close(fhand);
end