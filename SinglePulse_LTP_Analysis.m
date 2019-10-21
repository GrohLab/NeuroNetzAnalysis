%% LTP
% Loading the necessary files (spike times, 
% dataDir = 'E:\Data\VPM\LTP\190701_LTP_3700_1500_1520';
dataDir = 'D:\LTP\190716_Jesus_Emilio LTP_3751_1520_1500';
figureDir = fullfile(dataDir,'Figures\');
if ~loadTriggerData(dataDir)
    fprintf(1,'Not possible to load all the necessary variables')
    return
end
%% Constructing the helper 'global' variables
% Number of conditions
Nccond = numel(Conditions);
% Number of total samples
Ns = min(structfun(@numel,Triggers));
% Total duration of the recording
Nt = Ns/fs;
% Useless clusters (labeled as noise or they have very low firing rate)
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
bads = find(badsIdx);
totSpkCount = cellfun(@numel,sortedData(:,2));
clusterSpikeRate = totSpkCount/Nt;
silentUnits = clusterSpikeRate < 0.1;
bads = union(bads,find(silentUnits));
goods = setdiff(1:size(sortedData,1),bads);
badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));
% Logical spike trace for the first good cluster
spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
% Subscript column vectors for the rest good clusters
spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods(2:end),2),...
    'UniformOutput',false);
% Number of good clusters 
Ncl = numel(goods);
% Redefining the stimulus signals from the low amplitude to logical values
whStim = {'piezo','whisker'};
cxStim = {'laser','light'};
lfpRec = {'lfp','s1','cortex','s1lfp'};
trigNames = fieldnames(Triggers);
numTrigNames = numel(trigNames);
ctn = 1;
while ctn <= numTrigNames 
    if contains(trigNames{ctn},whStim,'IgnoreCase',true)
        whisker = Triggers.(trigNames{ctn})(1:Ns);
    end
    if contains(trigNames{ctn},cxStim,'IgnoreCase',true)
        laser = Triggers.(trigNames{ctn})(1:Ns);
    end
    if contains(trigNames{ctn},lfpRec,'IgnoreCase',true)
        LFP = Triggers.(trigNames{ctn})(1:Ns);
    end
    ctn = ctn + 1;
end
mObj = StepWaveform(whisker,fs);
mSubs = mObj.subTriggers;
piezo = mObj.subs2idx(mSubs,Ns);
lObj = StepWaveform(laser,fs);
lSubs = lObj.subTriggers;
laser = lObj.subs2idx(lSubs,Ns);
mObj.delete;lObj.delete;
continuousSignals = {double(piezo);double(laser);LFP};
%% User prompt for relevant information:
% Time lapse, bin size, and spontaneous and response windows
promptStrings = {'Viewing window (time lapse) [s]:','Response window [s]',...
    'Bin size [s]:'};
defInputs = {'0.1, 0.1', '0.002, 0.05', '0.01'};
answ = inputdlg(promptStrings,'Inputs', [1, 30],defInputs);
if isempty(answ)
    fprintf(1,'Cancelling...\n')
    return
else
    timeLapse = str2num(answ{1}); %#ok<*ST2NM>
    if numel(timeLapse) ~= 2
        timeLapse = str2num(inputdlg('Please provide the time window [s]:',...
            'Time window',[1, 30], '0.1, 0.1'));
        if isnan(timeLapse) || isempty(timeLapse)
            fprintf(1,'Cancelling...')
            return
        end
    end
    responseWindow = str2num(answ{2});
    binSz = str2double(answ(3));
end
fprintf(1,'Time window: %.2f - %.2f ms\n',timeLapse(1)*1e3, timeLapse(2)*1e3)
fprintf(1,'Response window: %.2f - %.2f ms\n',responseWindow(1)*1e3, responseWindow(2)*1e3)
fprintf(1,'Bin size: %.3f ms\n', binSz*1e3)
spontaneousWindow = -flip(responseWindow);
%% Condition triggered stacks
condNames = arrayfun(@(x) x.name,Conditions,'UniformOutput',false);
% Choose the conditions to create the stack upon
[chCond, iOk] = listdlg('ListString',condNames,'SelectionMode','single',...
    'PromptString','Choose the condition to look at:');
if ~iOk
    fprintf(1,'Cancelling...\n')
    return
end

% Select the onset or the offset of a trigger
fprintf(1,'Condition ''%s''\n', Conditions(chCond).name)
onOffStr = questdlg('Trigger on the onset or on the offset?','Onset/Offset',...
    'on','off','Cancel','on');
if strcmpi(onOffStr,'Cancel')
    fprintf(1,'Cancelling...\n')
    return
end

% Constructing the stack out of the user's choice
[dst, cst] = getStacks(spkLog,Conditions(chCond).Triggers,onOffStr,...
    timeLapse,fs,fs,spkSubs,continuousSignals);

[Ne, Nt, NTa] = size(dst);
% Computing the time axis for the stack
tx = (0:Nt)/fs - timeLapse(1);
%% LTP Exclusive: Control and post-induction conditions
% Boolean flags indicating which trigger belongs to which condition (delay
% flags)
% Determining if the gap between stimuli are big enough to separate the
% conditions
avGap = mean(diff(Conditions(chCond).Triggers(:,1),1,1));
[gapVal, timeGapSub] = max(diff(Conditions(chCond).Triggers(:,1),1,1));
errGap = log(avGap) - log(gapVal);
if abs(errGap) > 1
    % If the gap is bigger in one order of magnitude, we consider it to be
    % another condition (Control / Post-induction)
    condFlags = false(NTa, 2);
    condFlags(1:timeGapSub,1) = true;
    condFlags(timeGapSub+1:end,2) = true;
else
    % Otherwise, it is only one condition
    condFlags = true(NTa, 1);
end
Na = sum(condFlags,1);
%% Verification PSTH
consideredConditions = [1,3];
Nccond = 2;
goodsIdx = ~badsIdx';
for ccond = 1:Nccond
    [PSTH, trig, sweeps] = getPSTH(... dst([true;whiskerResponsiveUnitsIdx;true],:,:),timeLapse,...
        dst,timeLapse,...
        ~condFlags(:,ccond),binSz,fs);
    fig = plotClusterReactivity(PSTH,trig,sweeps,timeLapse,binSz,...
        [{Conditions(consideredConditions(ccond)).name};... sortedData(goods(whiskerResponsiveUnitsIdx),1);{'Laser'}],...
        sortedData(goods,1)],...
        strrep(expSubfix,'_','\_'));
    configureFigureToPDF(fig);
    print(fig,fullfile(figureDir,sprintf('%s %s.pdf',...
        expSubfix, Conditions(consideredConditions(ccond)).name)),...
        '-dpdf','-fillpage')
    print(fig,fullfile(figureDir,sprintf('%s %s.emf',...
        expSubfix, Conditions(consideredConditions(ccond)).name)),...
        '-dmeta')
end


%% Counting spikes in given windows and computing the statistical significance
sponActStackIdx = tx >= spontaneousWindow(1) & tx <= spontaneousWindow(2);
respActStackIdx = tx >= responseWindow(1) & tx <= responseWindow(2);
if size(condFlags,2) == 2
    % Getting the counts in the specified response and spontaneous time
    % windows
    contCountResponse = squeeze(sum(dst(2:end,respActStackIdx,condFlags(:,1)),2));
    contCountSpontan = squeeze(sum(dst(2:end,sponActStackIdx,condFlags(:,1)),2));
    postCountResponse = squeeze(sum(dst(2:end,respActStackIdx,condFlags(:,2)),2));
    postCountSpontan = squeeze(sum(dst(2:end,sponActStackIdx,condFlags(:,2)),2));
    % Computing the rates for the spike counts with the given time windows
    delta_t = diff(responseWindow);
    contRateResponse = mean(contCountResponse/delta_t,2);
    contRateSpontan = mean(contCountSpontan/delta_t,2);
    postRateResponse = mean(postCountResponse/delta_t,2);
    postRateSpontan = mean(postCountSpontan/delta_t,2);
    % Auxiliary variable to plot the diagonal line and to keep the scales
    % across the different plots 
    llx = max(max(max(contRateResponse, postRateResponse),...
        max(contRateSpontan, postRateSpontan)));
    % Computing the probability of spike for all trials in the given
    % windows
    % Flags
    contBoolResponse = double(contCountResponse > 0);
    contBoolSpontan = double(contCountSpontan > 0);
    postBoolResponse = double(postCountResponse > 0);
    postBoolSpontan = double(postCountSpontan > 0);
    % Probability
    contProbResponse = mean(contBoolResponse,2);
    contProbSpontan = mean(contBoolSpontan,2);
    postProbResponse = mean(postBoolResponse,2);
    postProbSpontan = mean(postBoolSpontan,2);
    
    % Shuffled counts for paired tests: 
    shufContCountResp = zeros(size(contCountResponse));
    shufContCountSpon = shufContCountResp;
    shufPostCountResp = zeros(size(postCountResponse));
    shufPostCountSpon = shufPostCountResp;
    for ccl = 1:Ncl 
        shufContCountResp(ccl,:) = contCountResponse(ccl,randperm(Na(1)));
        shufContCountSpon(ccl,:) = contCountSpontan(ccl,randperm(Na(1)));
        shufPostCountResp(ccl,:) = postCountResponse(ccl,randperm(Na(2)));
        shufPostCountSpon(ccl,:) = postCountSpontan(ccl,randperm(Na(2)));
    end
    
    % Computing the statistical tests: Wilcoxon [p, h] = signrank(x,y);
    % paired test, shuffling required!
    % Control and post-induction conditions for unaltered and shuffled
    % chonological sequence
    HContCond = zeros(Ncl,1);
    PContCond = HContCond;
    HShufContCond = HContCond;
    PShufContCond = HContCond;
    HPostCond = HContCond;
    PPostCond = HContCond;
    HShufPostCond = HContCond;
    PShufPostCond = HContCond;
    % KS test [h, p] = kstest2(x,y); for unpaired conditions 
    % Spontaneous and evoked activitiy in both control and post-induction
    % conditions.
    HSpon = HContCond;
    PSpon = HContCond;
    HResp = HContCond;
    PResp = HContCond;
    for ccl = 1:Ncl
        % Wilcoxon test
        % Control condition: spontaneous against evoked windows.
        [PContCond(ccl), HContCond(ccl)] = signrank(contCountSpontan(ccl,:),...
            contCountResponse(ccl,:));
        % Shuffled control condition
        [PShufContCond(ccl), HShufContCond(ccl)] = signrank(shufContCountSpon(ccl,:),...
            shufContCountResp(ccl,:));
        % Post-induction conditions
        [PPostCond(ccl), HPostCond(ccl)] = signrank(postCountSpontan(ccl,:),...
            postCountSpontan(ccl,:));
        % Shuffled post-induction condition
        [PShufPostCond(ccl), HShufPostCond(ccl)] = signrank(...
            shufPostCountSpon(ccl,:), shufPostCountResp(ccl,:));
        % Kolmogorov-Smirnov test
        % Spontaneous activity
        [HSpon(ccl), PSpon(ccl)] = kstest2(contCountSpontan(ccl,:),...
            postCountSpontan(ccl,:));
        % Evoked activity
        [HResp(ccl), PResp(ccl)] = kstest2(contCountResponse(ccl,:),...
            postCountResponse(ccl,:));
    end
else
end
%% Plotting the scatter data points
% Auxiliary variable for labelling the data points in the scatter plots.
clsLbls = num2str((1:Ncl)');
llx = round(llx*1.1,-1);
% Rate plots
controlFig = figure('Color',[1,1,1],'Name','Control condition','Visible','off');
postFig = figure('Color',[1,1,1],'Name','Post-induction','Visible','off');
sponFig = figure('Color',[1,1,1],'Name','Spontaneous','Visible','off');
respFig = figure('Color',[1,1,1],'Name','Evoked','Visible','off');
controlAx = axes('Parent',controlFig);
postAx = axes('Parent',postFig);
sponAx = axes('Parent',sponFig);
respAx = axes('Parent',respFig);

% Control plot (spontaneous vs evoked)
scatter(controlAx, contRateSpontan, contRateResponse, 'DisplayName',...
    'Rate per cluster'); grid(controlAx, 'on');
xlabel(controlAx, 'Spontaneous rate [Hz]');
ylabel(controlAx, 'Evoked rate [Hz]');
title(controlAx, 'Control condition: spontaneous vs evoked')
hold(controlAx,'on'); scatter(controlAx, contRateSpontan(HContCond == 1),...
    contRateResponse(HContCond == 1), 'Marker','square',...
    'DisplayName','Sync')
scatter(controlAx, contRateSpontan(HShufContCond == 1),...
    contRateResponse(HShufContCond == 1), 'Marker','+',...
    'DisplayName','Shuf')
legend(controlAx,'show'); axis(controlAx, 'square');  
line(controlAx, 'XData', [0, llx], 'YData', [0, llx], 'Color', [0.8, 0.8, 0.8],...
    'LineStyle', '--');
text(controlAx, contRateSpontan, contRateResponse, clsLbls);
controlFig.Visible = 'on';
% Post-induction (spontaneous vs evoked)
scatter(postAx, postRateSpontan, postRateResponse, 'DisplayName',...
    'Rate per cluster'); grid(postAx, 'on');
hold(postAx,'on'); scatter(postAx, postRateSpontan(HPostCond == 1),...
    postRateResponse(HPostCond == 1), 'Marker','square',...
    'DisplayName','Sync')
scatter(postAx, postRateSpontan(HShufPostCond == 1),...
    postRateResponse(HShufPostCond == 1), 'Marker','+',...
    'DisplayName','Shuf')
xlabel(postAx, 'Spontaneous rate [Hz]');
ylabel(postAx, 'Evoked rate [Hz]');
title(postAx, 'Post-induction condition: spontaneous vs evoked')
legend(postAx,'show');  axis(postAx, 'square');  
line(postAx, 'XData', [0, llx], 'YData', [0, llx], 'Color', [0.8, 0.8, 0.8],...
    'LineStyle', '--');
text(postAx, postRateSpontan, postRateResponse, clsLbls);
postFig.Visible = 'on';
% Spontaneous activity (Control vs post-induction)
scatter(sponAx, contRateSpontan, postRateSpontan); grid(sponAx, 'on');
xlabel(sponAx, 'Control spontaneous rate [Hz]');
ylabel(sponAx, 'Post-induction spontaneous rate [Hz]');
title(sponAx, 'Spontaneous activity: control vs post-induction')
hold(sponAx,'on'); scatter(sponAx, contRateSpontan(HSpon == 1),...
    postRateSpontan(HSpon == 1),'Marker','+','DisplayName','H=1');
legend(sponAx, 'show'); axis(sponAx, 'square');  
line(sponAx, 'XData', [0, llx], 'YData', [0, llx], 'Color', [0.8, 0.8, 0.8],...
    'LineStyle', '--');
text(sponAx, contRateSpontan, postRateSpontan, clsLbls);
sponFig.Visible = 'on';
% Evoked activity (Control vs post-induction)
scatter(respAx, contRateResponse, postRateResponse); grid(respAx, 'on');
xlabel(respAx, 'Control evoked rate [Hz]');
ylabel(respAx, 'Post-induction evoked rate [Hz]');
title(respAx, 'Evoked activity: control vs post-induction')
hold(respAx,'on'); scatter(respAx, contRateResponse(HResp == 1),...
    postRateResponse(HResp == 1),'Marker','+','DisplayName','H=1');
legend(respAx, 'show'); axis(respAx, 'square');  
line(respAx, 'XData', [0, llx], 'YData', [0, llx], 'Color', [0.8, 0.8, 0.8],...
    'LineStyle', '--');
text(respAx, contRateResponse, postRateResponse, clsLbls);
respFig.Visible = 'on';

% Probability plots
controlFigp = figure('Color',[1,1,1],'Name','Control condition','Visible','off');
postFigp = figure('Color',[1,1,1],'Name','Post-induction','Visible','off');
sponFigp = figure('Color',[1,1,1],'Name','Spontaneous','Visible','off');
respFigp = figure('Color',[1,1,1],'Name','Evoked','Visible','off');
controlAxp = axes('Parent',controlFigp);
postAxp = axes('Parent',postFigp);
sponAxp = axes('Parent',sponFigp);
respAxp = axes('Parent',respFigp);


% Control plot (spontaneous vs evoked)
scatter(controlAxp, contProbSpontan, contProbResponse, 'DisplayName',...
    'p(s) per cluster'); grid(controlAxp, 'on');
xlabel(controlAxp, 'p(Spontaneous)');
ylabel(controlAxp, 'p(Evoked)');
title(controlAxp, 'Control condition_{p(s)}: spontaneous vs evoked')
%{ 
% STATISTICAL SIGNIFICANCE HIGHLIGHTING
% hold(controlAxp,'on'); scatter(controlAxp, contProbSpontan(HContCond == 1),...
%     contProbResponse(HContCond == 1), 'Marker','square',...
%     'DisplayName','Sync')
% scatter(controlAxp, contProbSpontan(HShufContCond == 1),...
%     contProbResponse(HShufContCond == 1), 'Marker','+',...
%     'DisplayName','Shuf')
% legend(controlAxp,'show')
%}
text(controlAxp, contProbSpontan, contProbResponse, clsLbls);
axis(controlAxp, 'square');  
line(controlAxp, 'XData', [0, 1], 'YData', [0, 1], 'Color', [0.8, 0.8, 0.8],...
    'LineStyle', '--');
controlFigp.Visible = 'on';
% Post-induction (spontaneous vs evoked)
scatter(postAxp, postProbSpontan, postProbResponse, 'DisplayName',...
    'p(s) per cluster'); grid(postAxp, 'on');
%{
STATISTICAL SIGNIFICANCE HIGHLIGHTING
hold(postAxp,'on'); scatter(postAxp, postProbSpontan(HPostCond == 1),...
    postProbResponse(HPostCond == 1), 'Marker','square',...
    'DisplayName','Sync')
scatter(postAxp, postProbSpontan(HShufPostCond == 1),...
    postProbResponse(HShufPostCond == 1), 'Marker','+',...
    'DisplayName','Shuf')
legend(postAxp,'show')
%}
xlabel(postAxp, 'p(Spontaneous)');
ylabel(postAxp, 'p(Evoked)');
title(postAxp, 'Post-induction condition_{p(s)}: spontaneous vs evoked')
text(postAxp, postProbSpontan, postProbResponse, clsLbls);
axis(postAxp, 'square');  
line(postAxp, 'XData', [0, 1], 'YData', [0, 1], 'Color', [0.8, 0.8, 0.8],...
    'LineStyle', '--');
postFigp.Visible = 'on';
% Spontaneous activity (Control vs post-induction)
scatter(sponAxp, contProbSpontan, postProbSpontan); grid(sponAxp, 'on');
xlabel(sponAxp, 'p(Spontaneous | Control)');
ylabel(sponAxp, 'p(Spontaneous | Post-induction)');
title(sponAxp, 'Spontaneous activity_{p(s)}: control vs post-induction')
%{
STATISTICAL SIGNIFICANCE HIGHLIGHTING
hold(sponAxp,'on'); scatter(sponAxp, contProbSpontan(HSpon == 1),...
    postProbSpontan(HSpon == 1),'Marker','+','DisplayName','H=1');
%}    
text(sponAxp, contProbSpontan, postProbSpontan, clsLbls);
axis(sponAxp, 'square');  
line(sponAxp, 'XData', [0, 1], 'YData', [0, 1], 'Color', [0.8, 0.8, 0.8],...
    'LineStyle', '--');
sponFigp.Visible = 'on';
% Evoked activity (Control vs post-induction)
scatter(respAxp, contProbResponse, postProbResponse); grid(respAxp, 'on');
xlabel(respAxp, 'p(Evoked | Control)');
ylabel(respAxp, 'p(Evoked | Post-induction)');
title(respAxp, 'Evoked activity_{p(s)}: control vs post-induction')
%{
STATISTICAL SIGNIFICANCE HIGHLIGHTING
hold(respAxp,'on'); scatter(respAxp, contProbResponse(HResp == 1),...
    postProbResponse(HResp == 1),'Marker','+','DisplayName','H=1');
%}
text(respAxp, contProbResponse, postProbResponse, clsLbls);
axis(respAxp, 'square');  
line(respAxp, 'XData', [0, 1], 'YData', [0, 1], 'Color', [0.8, 0.8, 0.8],...
    'LineStyle', '--');
respFigp.Visible = 'on';
%% Saving the results (Figures and configuration)
pause(20);
saveFigs =...
    questdlg('Do you wish to save the figures?','Save','Yes','No','Yes');
saveFlag = false;
if strcmpi(saveFigs, 'yes')
    saveFlag = true;
    if ~exist(figureDir,'dir')
        if ~mkdir(figureDir)
            fprintf(1,'There was a problem creating the ''Figures'' folder\n');
            fprintf(1,'Please verify the location\n')
            saveFlag = false;
        end
    end
else
    closeFigs = questdlg('Would you like to close all figures?',...
        'Close all', 'Yes', 'No', 'No');
    if strcmpi(closeFigs,'yes')
        close all
    end
end

saveConfig = questdlg('What about the configuration for the analysis?',...
    'Save Configuration', 'Yes', 'No', 'Yes');
if strcmpi(saveConfig,'yes')
    configStruct = struct('TimeLapse',timeLapse,...
        'ResponseWindow',responseWindow,...
        'BinSize',binSz, 'Condition', Conditions(chCond).name);
    configFile = [expSubfix,'_ConfigStruct.mat'];
    if ~exist(configFile,'file')
        save([expSubfix,'_ConfigStruct.mat'],'configStruct')
    else
        fprintf(1,'The configuration file exists! No configuration saved\n')
    end
end

if saveFlag
    figArray = [controlFig,postFig,sponFig,respFig,...
        controlFigp,postFigp,sponFigp,respFigp];
    figNameCell = {'Control condition', 'Post-induction condition',...
        'Spontaneous activity', 'Evoked activity'};
    figNameCell = [figNameCell, ...
        cellfun(@(x) [x,' p(s)'],figNameCell,'UniformOutput',false)];
    arrayfun(@configureFigureToPDF, figArray);
    for cf = 1:numel(figArray)
        emfFile = fullfile(figureDir,[figNameCell{cf},'.emf']);
        pdfFile = fullfile(figureDir,[figNameCell{cf},'.pdf']);
        figFile = fullfile(figureDir,[figNameCell{cf},'.fig']);
        fileExist = [exist(emfFile,'file'),...
            exist(pdfFile,'file'), exist(figFile,'file')];
        if any(fileExist)
            ovrAns =...
                questdlg('Some figures exist already. Do you wish to overwrite?',...
                'Yes', 'No', 'No');
            if strcmp(ovrAns,'No')                            
                fprintf(1,'No figures saved!\n')
                break;
            end               
        end
        print(figArray(cf), emfFile, '-dmeta')
        print(figArray(cf), pdfFile, '-dpdf', '-fillpage')
        savefig(figArray(cf), figFile)
    end
end
