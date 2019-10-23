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

binFiles = dir(fullfile(dataDir,'*.bin'));
[~, baseName, ~] = fileparts(binFiles(1).name);

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
    existFlags = false;
    for ch = 1:numel(headFields)
        if contains(headVars.(headFields{ch}).title, whiskStimClue,...
                'IgnoreCase', true)
            existFlags = true;
            wChanField = strrep(headFields{ch},'head','chan');
            currWhisk = chanVars.(wChanField);
            if iscolumn(currWhisk)
                currWhisk = currWhisk';
            end
            whisk = [whisk, currWhisk(1:N)];
        end
        if contains(headVars.(headFields{ch}).title, laserStimClue,...
                'IgnoreCase', true)
            existFlags = false;
            lChanField = strrep(headFields{ch},'head','chan');
            currLaser = chanVars.(lChanField);
            if iscolumn(currLaser)
                currLaser = currLaser';
            end
            laser = [laser, currLaser(1:N)];
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
    if existFlags
        if isempty(laser)
            laser = [laser, repmat(min(whisk),1,N)];
        else
            laser = [laser, repmat(min(laser),1,N)];
        end
    else
        whisk = [whisk, repmat(min(whisk),1,N)];
    end
end

%% Computation of the subscripts and creation of the condition variable
pObj = StepWaveform(whisk, fs);
pSubs = pObj.subTriggers;
if any(diff(pSubs(:,1),1,1) < 1)
    pFrstFlag = pObj.FirstOfTrain;
    pSubs = pSubs(pFrstFlag,:);
end
% Big silence in between control and post-induction
[~, condSepSub] = max(diff(pSubs(:,1),1,1));


lObj = StepWaveform(laser, fs);
lSubs = lObj.subTriggers;
lFrstFlag = lObj.FirstOfTrain;
lSubs = lSubs(lFrstFlag,:);

Conditions = struct('name', 'Control',...
    'Triggers', pSubs(1:condSepSub,:));
Conditions(2).name = 'Induction';
Conditions(2).Triggers = lSubs;
Conditions(3).name = 'Post_induction';
Conditions(3).Triggers =  pSubs(condSepSub+1:end,:);
Conditions(4).name = 'AllPiezo';
Conditions(4).Triggers = pSubs;
Triggers = struct('Piezo', whisk, 'Laser', laser, 'LFP', lfp);

analysisFileName = fullfile(dataDir,[baseName, 'analysis.mat']);
if exist(analysisFileName,'file')
    ovwrAns = questdlg(['The analysis file for the given folder exists.',...
        ' Would you like to overwrite it?'],'Overwrite','Yes','No','Yes');
    if strcmp(ovwrAns,'No')
        fprintf(1,'No file saved.')
        return
    end
end
save(analysisFileName, 'Conditions', 'Triggers', 'fs');

end
