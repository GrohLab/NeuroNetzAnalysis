function mBehRes = analyseBehaviour_allSessions(expDirs)

fnOpts = {'UniformOutput', false};
getChildFolder = @(x) fullfile(x.folder, x.name);

mBehRes = cell(numel(expDirs),1);
delSub = [];
for ec = 1:numel(expDirs)
    % Load *analysis.mat file from the ephys_* folder.
    cndFile = dir(fullfile(getChildFolder(expDirs(ec)), "*\*analysis.mat"));
    cmts = dir(fullfile(getChildFolder(expDirs(ec)), 'Comments.txt'));
    if isempty(cndFile) || ~isempty(cmts)
        fprintf(1, 'Experiment not analysed!\n')
        [~, mouseName] = fileparts(expDirs(ec).folder);
        fprintf(1, '%s\n', [mouseName, ' ', expDirs(ec).name])
        delSub = [delSub, ec];
        continue
    end
    load(getChildFolder(cndFile), 'Conditions')

    pairedStim = arrayfun(@(x) Conditions(1).Triggers(:,1) == ...
        Conditions(x).Triggers(:,1)', 2:numel(Conditions), fnOpts{:});
    pairedStim = cellfun(@(x) any(x, 2), pairedStim, fnOpts{:});
    pairedStim = cat(2, pairedStim{:});
    consCondNames = arrayfun(@(c) string(Conditions(c).name), 1+find(sum(pairedStim)));
    pairedStim(:,sum(pairedStim)==0) = [];

    behDir = fullfile(getChildFolder(expDirs(ec)), 'Behaviour');
    figureDir = recursiveFolderSearch(getChildFolder(expDirs(ec)), 'Figures');
    if ~isempty(figureDir)
        [behRes, behFigDir] = analyseBehaviour(behDir, 'Condition', 'P', ...
            'PairedFlags', pairedStim, 'FigureDirectory', figureDir, ...
            'ConditionsNames', cellstr(consCondNames), 'verbose', false, ...
            'showPlots', false);

        biFigPttrn = "BehIndex%s";
        biFigPttrn = sprintf(biFigPttrn, sprintf(" %s (%%.3f)", consCondNames));
        [pAreas, ~, behAreaFig] = createBehaviourIndex(behRes);
        behRes = arrayfun(@(bs, ba) setfield(bs,'BehIndex', ba), behRes, pAreas);
        set(behAreaFig, 'UserData', behRes)

        biFN = sprintf(biFigPttrn, pAreas);

        saveFigure(behAreaFig, fullfile(behFigDir, biFN), true, true);
        close all;
        mBehRes{ec} = behRes;
    end
end
mBehRes(delSub) = [];