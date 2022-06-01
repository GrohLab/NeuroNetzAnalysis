function [outputArg1,outputArg2] = analyseBehaviour(behDir, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Auxiliary variables
afPttrn = "ArduinoTriggers*.mat";
rfPttrn = "RollerSpeed*.mat";
lPttrn = "Laser*.csv"; pPttrn = "Puff*.csv";
fnOpts = {'UniformOutput', false};
axOpts = {'Box','off','Color','none'};
lgOpts = cat(2, axOpts{1:2}, {'Location','best'});
m = 1e-3; k = 1/m;
%% Input validation
p = inputParser;

defVW = [-250, 500]*m; %milliseconds
checkVW = @(x) all([isnumeric(x), numel(x) == 2, x(1)<x(2)]);

defRW = [5, 400]*m;
checkTH = @(x) x > 0 & numel(x) == 1;
defSig = 2.5; % Spontaneous standard deviation
defSMed = 0.2; % Spontaneous Median
defMed = 1; % Median in all viewing window

defSpeedTh = 0.1:0.1:3; % Speed thresholds
checkSpeedTh = @(x) all([diff(x) > 0, x > 0, numel(x) > 1]);

addRequired(p, 'behDir', @(x) ischar(x) | isstring(x))
addParameter(p, 'ViewingWindow', defVW, checkVW)
addParameter(p, 'ResponseWindow', defRW, checkVW)
addParameter(p, 'SpeedTh', defSpeedTh, checkSpeedTh)
addParameter(p, 'SponSigTh', defSig, checkTH)
addParameter(p, 'SponMedianTh', defSMed, checkTH)
addParameter(p, 'MedianTh', defMed, checkTH)

parse(p, behDir, varargin{:})

behDir = p.Results.behDir;
bvWin = p.Results.ViewingWindow;
brWin = p.Results.ResponseWindow;
spTh = {p.Results.SpeedTh}; 
sigTh = p.Results.SponSigTh; 
sMedTh = p.Results.SponMedianTh;
tMedTh = p.Results.MedianTh;

if bvWin(1) > brWin(1) || bvWin(2) < brWin(2)
    fprintf(1, 'Given response window lays outside the viewing window!\n')
    fprintf(1, 'VW:%.2f - %.2f ms\nRW:%.2f - %.2f ms\n', bvWin*k, brWin*k)
    fprintf(1, 'Resetting to default windows!\n')
    bvWin = defVW; brWin = defRW;
    fprintf(1, 'VW:%.2f - %.2f ms\nRW:%.2f - %.2f ms\n', bvWin*k, brWin*k)
end

if ~exist(behDir, "dir")
    fprintf(1, 'Given folder doesn''t exist!\nExiting...')
    return
end

%%
% If only one folder named Behaviour exists, chances are that this is
% an awake experiment.
if isempty(dir(fullfile(behDir, afPttrn)))
    readAndCorrectArdTrigs(behDir);
end
fprintf(1,'Time window: %.2f - %.2f ms\n',bvWin*k)
fprintf(1,'Response window: %.2f - %.2f ms\n',brWin*k)
% Roller speed
rfFiles = dir(fullfile(behDir, rfPttrn));
if isempty(rfFiles)
    [~, vf, ~, fr, Texp] = createRollerSpeed(behDir);
    rfFiles = dir(fullfile(behDir, rfPttrn));
end
if numel(rfFiles) == 1
    rfName = fullfile(rfFiles.folder, rfFiles.name);
    load(rfName) %#ok<LOAD>
    try
        %              Encoder steps  Radius^2
        en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*rollFs;
        fr = rollFs;
    catch
        try
            %              Encoder steps  Radius^2
            en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*fr;
        catch
            en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*fsRoll;
            fr = fsRoll;
        end
    end
end
% Triggers
getFilePath = @(x) fullfile(x.folder, x.name);
atVar = {'atTimes', 'atNames', 'itTimes', 'itNames'};
afFiles = dir(fullfile(behDir,afPttrn));
if ~isempty(afFiles)
    atV = arrayfun(@(x) load(getFilePath(x), atVar{:}), afFiles);
    atT = arrayfun(@(x, z) cellfun(@(y, a) y+a, ...
        x.atTimes, repmat(z,1,length(x.atTimes)), fnOpts{:}), atV', ...
        num2cell([0, Texp(1:end-1)]), fnOpts{:});
    atT = cat(1, atT{:});
    atTimes = arrayfun(@(x) cat(1, atT{:,x}), 1:size(atT,2), fnOpts{:});
    atNames = atV(1).atNames;
end

lSub = arrayfun(@(x) contains(Conditions(chCond).name, x), atNames);
[~, vStack] = getStacks(false, round(atTimes{lSub} * fr), 'on', bvWin,...
    fr, fr, [], vf*en2cm); [Nbt, Ntr] = size(vStack, [2, 3]); Na = Ntr; delayFlags = true(Na, 1);
tmdl = fit_poly([1,Nbt], bvWin, 1);
behTx = ((1:Nbt)'.^[1,0])*tmdl;
% Spontaneous flag
bsFlag = behTx <= 0; brFlag = behTx < brWin;
brFlag = xor(brFlag(:,1),brFlag(:,2));
sSig = squeeze(std(vStack(:,bsFlag,:), [], 2));
sMed = squeeze(median(vStack(:,bsFlag,:), 2));
tMed = squeeze(median(vStack, 2));

% A bit arbitrary threshold, but enough to remove running trials
thrshStr = sprintf("TH s%.2f sp_m%.2f t_m%.2f", sigTh, sMedTh, tMedTh);
excFlag = sSig > sigTh | abs(sMed) > sMedTh | abs(tMed) > tMedTh;
ptOpts = {"Color", 0.7*ones(1,3), "LineWidth", 0.2;...
    "Color", "k", "LineWidth",  1.5};

gp = zeros(Nccond, 1, 'single');
rsPttrn = "%s roller speed VW%.2f - %.2f s RM%.2f - %.2f ms EX%d %s";
pfPttrn = "%s move probability %.2f RW%.2f - %.2f ms EX%d %s";
rsSgnls = cell(Nccond, 1); mvFlags = cell(Nccond,1); mvpt = mvFlags;
qSgnls = rsSgnls;
mat2ptch = @(x) [x(1:end,:)*[1;1]; x(end:-1:1,:)*[1;-1]];
getThreshCross = @(x) sum(x)/size(x,1);
xdf = arrayfun(@(x) ~excFlag & delayFlags(:,x), 1:Nccond, ...
    fnOpts{:});  xdf = cat(2, xdf{:});

for ccond = 1:Nccond
    sIdx = xdf(:,ccond);
    % % Plot speed signals
    fig = figure("Color", "w");
    Nex = sum(xor(sIdx, delayFlags(:,ccond)));
    rsFigName = sprintf(rsPttrn,consCondNames{ccond}, bvWin,...
        brWin*1e3, Nex, thrshStr);
    % Plot all trials
    plot(behTx, squeeze(vStack(:,:,sIdx)), ptOpts{1,:}); hold on;
    % Plot mean of trials
    % Standard deviation
    %rsSgnls{ccond} = [squeeze(mean(vStack(:,:,sIdx),3))',...
    %squeeze(std(vStack(:,:,sIdx),1,3))'];
    % S.E.M.
    rsSgnls{ccond} = [squeeze(mean(vStack(:,:,sIdx),3))',...
        squeeze(std(vStack(:,:,sIdx),1,3))'./sqrt(sum(sIdx))];
    qSgnls{ccond} = squeeze(quantile(vStack(:,:,sIdx),3,3));
    lObj = plot(behTx, rsSgnls{ccond}(:,1), ptOpts{2,:});
    lgnd = legend(lObj,string(consCondNames{ccond}));
    set(lgnd, "Box", "off", "Location", "best")
    set(gca, axOpts{:})
    title(['Roller speed ',consCondNames{ccond}])
    xlabel("Time [s]"); ylabel("Roller speed [cm/s]"); xlim(bvWin)
    saveFigure(fig, fullfile(figureDir, rsFigName), 1)
    % Probability plots
    mvpt{ccond} = getMaxAbsPerTrial(squeeze(vStack(:,:,sIdx)), ...
        brWin, behTx);
    mvFlags{ccond} = compareMaxWithThresh(mvpt{ccond}, spTh);
    gp(ccond) = getAUC(mvFlags{ccond});
    pfName = sprintf(pfPttrn, consCondNames{ccond}, gp(ccond),...
        brWin*1e3, Nex, thrshStr);
    fig = plotThetaProgress(mvFlags(ccond), spTh,...
        string(consCondNames{ccond}));
    xlabel("Roller speed \theta [cm/s]");
    title(sprintf("Trial proportion crossing \\theta: %.3f", gp(ccond)))
    saveFigure(fig, fullfile(figureDir, pfName), 1)
end
clMap = lines(Nccond);
phOpts = {'EdgeColor', 'none', 'FaceAlpha', 0.25, 'FaceColor'};
% Plotting mean speed signals together
fig = figure("Color", "w"); axs = axes("Parent", fig, "NextPlot", "add");
arrayfun(@(x) patch(axs, behTx([1:end, end:-1:1]),...
    mat2ptch(rsSgnls{x}), 1, phOpts{:}, clMap(x,:)), 1:Nccond); hold on
lObj = arrayfun(@(x) plot(axs, behTx, rsSgnls{x}(:,1), "Color", clMap(x,:),...
    "LineWidth", 1.5, "DisplayName", consCondNames{x}), 1:Nccond);
xlabel(axs, "Time [s]"); xlim(axs, bvWin); ylabel(axs, "Roller speed [cm/s]")
set(axs, axOpts{:}); title(axs, "Roller speed for all conditions")
lgnd = legend(axs, lObj); set(lgnd, lgOpts{:})
rsPttrn = "Mean roller speed %s VW%.2f - %.2f s RM%.2f - %.2f ms EX%s %s SEM";
Nex = Na - sum(xdf);
rsFigName = sprintf(rsPttrn, sprintf('%s ', consCondNames{:}), bvWin,...
    brWin*1e3, sprintf('%d ', Nex), thrshStr);
saveFigure(fig, fullfile(figureDir, rsFigName), 1)

% Plotting median speed signals together
q2patch = @(x) [x(:,1);x(end:-1:1,3)];
fig = figure("Color", "w"); axs = axes("Parent", fig, "NextPlot", "add");
arrayfun(@(x) patch(axs, behTx([1:end, end:-1:1]),...
    q2patch(qSgnls{x}), 1, phOpts{:}, clMap(x,:)), 1:Nccond); hold on
lObj = arrayfun(@(x) plot(axs, behTx, qSgnls{x}(:,2), "Color", clMap(x,:),...
    "LineWidth", 1.5, "DisplayName", consCondNames{x}), 1:Nccond);
xlabel(axs, "Time [s]"); xlim(axs, bvWin); ylabel(axs, "Roller speed [cm/s]")
set(axs, axOpts{:}); title(axs, "Roller speed for all conditions")
lgnd = legend(axs, lObj); set(lgnd, lgOpts{:})
rsPttrn = "Median roller speed %s VW%.2f - %.2f s RM%.2f - %.2f ms EX%s %s IQR";
rsFigName = sprintf(rsPttrn, sprintf('%s ', consCondNames{:}), bvWin,...
    brWin*1e3, sprintf('%d ', Nex), thrshStr);
saveFigure(fig, fullfile(figureDir, rsFigName), 1)

% Plotting movement threshold crossings
fig = figure("Color", "w"); axs = axes("Parent", fig, "NextPlot", "add");
mvSgnls = cellfun(getThreshCross, mvFlags, fnOpts{:});
mvSgnls = cat(1, mvSgnls{:}); mvSgnls = mvSgnls';
plot(axs, spTh{1}, mvSgnls);
ccnGP = cellfun(@(x, y) [x, sprintf(' AUC%.3f',y)], consCondNames', ...
    num2cell(gp), fnOpts{:});
lgnd = legend(axs, ccnGP); set(axs, axOpts{:})
set(lgnd, lgOpts{:}); ylim(axs, [0,1])
xlabel(axs, "Roller speed \theta [cm/s]"); ylabel(axs, "Trial proportion")
title(axs, "Trial proportion crossing \theta")
pfPttrn = "Move probability %sRW%.2f - %.2f ms %s";
pfName = sprintf(pfPttrn, sprintf('%s ', ccnGP{:}), brWin*1e3, thrshStr);
saveFigure(fig, fullfile(figureDir, pfName), 1)
% Tests for movement
prms = nchoosek(1:Nccond,2);
getDistTravel = @(x) squeeze(sum(abs(vStack(:,brFlag,xdf(:,x))),2));
dstTrav = arrayfun(getDistTravel, 1:Nccond, fnOpts{:});
[pd, hd, statsd] = arrayfun(@(x) ranksum(dstTrav{prms(x,1)}, ...
    dstTrav{prms(x,2)}), 1:size(prms,1), fnOpts{:});
[pm, hm, statsm] = arrayfun(@(x) ranksum(mvpt{prms(x,1)}, ...
    mvpt{prms(x,2)}), 1:size(prms,1), fnOpts{:});
%%
end