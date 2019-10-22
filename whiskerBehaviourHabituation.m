%% Whisker adaptation behaviour
% Running through the whisker files of the wild type animals and observing
% their whisking patterns through time. We are looking for whisking
% behaviour habituation.
m = 1e-3; k = 1e3;
fs1 = k;
fs20 = 20*k;
frstCell = @(x) x{1}; frstRow = @(x) x(1,:);
%% Directory selection
% Directories selection containing the whisker's information (whiskDir) and
% the triggered time points (analysisDir)
fprintf(1,'Please select a WHISKER information directory\n')
whiskDir = uigetdir('Whisker directory');
if ~ischar(whiskDir)
    fprintf(1,'Cancelling...\n')
    return
end

fprintf(1,'Please select a TRIGGER directory\n')
analysisDir = uigetdir(getParentDir(whiskDir,1),'Analysis directory');
if ~ischar(analysisDir)
    fprintf(1,'Cancelling...\n')
    return
end

% Whisking and analysis files
whiskFiles = dir(fullfile(whiskDir,'T*.mat'));
analysisFiles = dir(fullfile(analysisDir,'*\T*.mat'));

% Number of files found in both directories
Nwf = numel(whiskFiles);
Naf = numel(analysisFiles);
% Verification if the file numbers are the same
if Nwf ~= Naf
    fprintf(1,'Please verify the files in the given folders!')
    return
end

% Determining the animal IDs
animalIDsCell = ...
    cellfun(frstCell,(cellfun(@strsplit,...
    frstRow(struct2cell(whiskFiles))', repmat({'_'},Nwf,1),...
    'UniformOutput',false)),'UniformOutput',false);
animalIDs = unique(cellfun(@str2double, cellfun(@strrep, animalIDsCell,...
    repmat({'T'},Nwf,1),repmat({''},Nwf,1),'UniformOutput', false)));

%% Main loop
timeLapse = [0.5, 5];
figureDir = fullfile(getParentDir(whiskDir,2),'Figures Whisker');
for canimal = 1:numel(animalIDs)
    cst = [];
    fprintf(1,'Animal %d ',animalIDs(canimal))
    animalFilesFlag = cellfun(@contains, animalIDsCell,...
        repmat({num2str(animalIDs(canimal))},Nwf,1));
    fprintf(1,'has %d file(s)\n',sum(animalFilesFlag))
    files2process = find(animalFilesFlag)';
    Na = zeros(sum(animalFilesFlag),1,'int16');
    caux = 1;
    for cfs = files2process
        load(fullfile(whiskFiles(cfs).folder, whiskFiles(cfs).name),...
            'WhiskerAngle')
        load(fullfile(analysisFiles(cfs).folder, analysisFiles(cfs).name),...
            'Triggers')
        pObj = StepWaveform(double(Triggers.puff), fs1);
        pSubs = pObj.subTriggers;
        Ns = pObj.NSamples;
        % Computing a semi-periodic pulse flag.
        puffDist = distmatrix(pSubs(:,1)/fs20,pSubs(:,1)/fs20);
        Nta = size(puffDist,1);
        cp = 1;
        periodicPulseFlag = [true;false(Nta-1,1)];
        while cp <= Nta
            cp = find(puffDist(cp:end,cp) >= 5,1,'first') + cp - 1;
            periodicPulseFlag(cp) = true;
        end
        fprintf(1,'Periodic pulses found in %s: %d\n',whiskFiles(cfs).name,...
            sum(periodicPulseFlag))
        pSubs = pSubs(periodicPulseFlag,:);
        Na(caux) = size(pSubs,1);
        [~, cstaux] = getStacks(false(1,Ns), pSubs, 'on', timeLapse, fs20,...
            fs1, [], {WhiskerAngle});
        cst = cat(3, cst, cstaux);
        caux = caux + 1;
    end
    whiskBehave = squeeze(cst)';
    % Figure creation and saving
    whiskBehaveNorm = whiskBehave./max(whiskBehave,[],2);
    tx1st = (0:size(cst,2)-1) * 1/fs1 - timeLapse(1);
    fig = figure('Name', sprintf('Whisking behaviour for animal %s',...
        num2str(animalIDs(canimal))), 'Visible', 'off');
    ax(1) = subplot(6,1,1:5,'Parent',fig);
    imagesc(ax(1), 'XData', tx1st, 'YData', 1:sum(Na),...
        'CData', whiskBehaveNorm)
    ax(1).XAxis.Visible = 'off';
    ax(2) = subplot(6,1,6, 'Parent', fig);
    avWhiskBeh = mean(whiskBehave,1);
    plot(ax(2), tx1st, avWhiskBeh, 'DisplayName', 'Mean whisker angle')
    legend(ax(2), 'show')
    linkaxes(ax, 'x')
    axis(ax(2), [-timeLapse(1), timeLapse(2),...
        min(avWhiskBeh)*0.99, max(avWhiskBeh)*1.01])
    axis(ax(1), [-timeLapse(1), timeLapse(2), 0.5, sum(Na) + 0.5])
    title(ax(1), ['Whisking behaviour for animal ',...
        num2str(animalIDs(canimal))])
    xlabel(ax(2), 'Time [s]')
    ylabel(ax(1), 'Trial number_{Experiment division}')
    ylabel(ax(2), '\angle Whisker [?]')
    fig.Visible = 'on';
    configureFigureToPDF(fig);
    pause(5)
    saveAns = questdlg(['Save figure for animal ',...
        num2str(animalIDs(canimal))], 'Save figure' ,'Yes', 'No', 'Yes');
    if strcmp(saveAns,'Yes')
        figureFileName = fullfile(figureDir,...
            sprintf('Whisking behaviour T%d',animalIDs(canimal)));
        savefig(fig, [figureFileName,'.fig'])
        print(fig, [figureFileName, '.emf'], '-dmeta')
        print(fig, [figureFileName, '.pdf'], '-dpdf', '-fillpage')
    end
end