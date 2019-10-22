function mergeTriggers(dataDir)
% MERGETRIGGERS collects the triggers signals of several VPM experiments
% and merges them into one analysis file.

frstRow = @(x) x(1,:);
confSigFiles = dir(fullfile(dataDir,'*_CondSig.mat'));
if isempty(confSigFiles)
    confSigFiles = dir(fullfile(dataDir,'*\*_CondSig.mat'));
end
if isempty(confSigFiles)
    fprintf(1,'Didn''t find the necessary files in this folder.\n')
    fprintf(1,'Please verify the existance of *_CondSig.mat file\n')
    return
end
Ncf = numel(confSigFiles);
% Selecting files
fileNamesCell = frstRow(struct2cell(confSigFiles));
[incFiles, iok] = listdlg('ListString',fileNamesCell,...
    'SelectionMode','multiple',...
    'PromptString','Select the files to join:',...
    'InitialValue',1:Ncf);
if iok
    confSigFiles = confSigFiles(incFiles);
    fileNamesCell = frstRow(struct2cell(confSigFiles));
else
    fprintf(1,'Cancelling...\n')
    return
end
Ncf = numel(confSigFiles);
% Getting the correct order
fileOrder = (1:Ncf)';
defInput = num2cell(num2str(fileOrder));
answr = inputdlg(fileNamesCell, 'File order', [1, 60], defInput);
if isempty(answr)
    fprintf(1,'Cancelling...\n')
    return
end
nFileOrder = str2double(answr);
nConfSigFiles = confSigFiles;
if ~isempty(answr) && sum(abs(fileOrder - nFileOrder)) ~= 0
    fprintf(1,'Changing file order...\n')
    nConfSigFiles(nFileOrder) = confSigFiles(fileOrder);
    confSigFiles = nConfSigFiles;
else
    fprintf('File order not altered\n')
end
%% Main loop

whiskStimClue = {'piezo', 'whisker', 'puff'};
laserStimClue = {'light', 'laser'};
cortxRecClue = {'lfp', 'ctx'};
whisk = []; laser = []; lfp = [];
for cfs = 1:Ncf
    fprintf(1,'%s\n', confSigFiles(cfs).name)
    auxFile = fullfile(confSigFiles(cfs).folder, confSigFiles(cfs).name);
    headVars = load(auxFile, 'head*');
    chanVars = load(auxFile, 'chan*');
    headFields = fieldnames(headVars);
    N = max(structfun(@(x) x.npoints,headVars));
    fs = max(structfun(@(x) x.SamplingFrequency,headVars));
    existFlags = false(1,2);
    for ch = 1:numel(headFields)
        if contains(headVars.(headFields{ch}).title, whiskStimClue,...
                'IgnoreCase', true)
            existFlags(1) = true;
            wChanField = strrep(headFields{ch},'head','chan');
            currWhisk = chanVars.(wChanField);
            if iscolumn(currWhisk)
                currWhisk = currWhisk';
            end
            whisk = [whisk, currWhisk];
        end
        if contains(headVars.(headFields{ch}).title, laserStimClue,...
                'IgnoreCase', true)
            existFlags(2) = true;
            lChanField = strrep(headFields{ch},'head','chan');
            currLaser = chanVars.(lChanField);
            if iscolumn(currLaser)
                currLaser = currLaser';
            end
            laser = [laser, currLaser];
        end
        if contains(headVars.(headFields{ch}).title, cortxRecClue,...
                'IgnoreCase', true)
            cChanField = strrep(headFields{ch},'head','chan');
            lfs = headVars.(headFields{ch}).SamplingFrequency;
            [currlfp, tx] = changeSignalFS(chanVars.(cChanField), lfs, fs);
            if iscolumn(currlfp)
                currlfp = currlfp';
            end
            try
                lfp = [lfp, currlfp(1:N)];
            catch
                auxlfp = padarray(currlfp,[0,double(N-length(currlfp))],...
                    'symmetric','post');
                lfp = [lfp, auxlfp];
            end
        end
    end
    if existFlags(1)
        if isempty(laser)
            laser = [laser, repmat(min(whisk),1,N)];
        else
            laser = [laser, repmat(min(laser),1,N)];
        end
    else
        whisk = [whisk, repmat(min(whisk),1,N)];
    end
end
fprintf(1,'To be continued...\n')

end
