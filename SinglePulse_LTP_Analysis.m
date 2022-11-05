%% LTP
% Choosing the working directory
dataDir = uigetdir('E:\Data\VPM\LTP','Choose a working directory');
if dataDir == 0
    return
end
%dataDir = 'E:\Data\VPM\LTP\191016_Jesus_LTP_3710_1520_1500';
% dataDir = 'D:\LTP\190716_Jesus_Emilio LTP_3751_1520_1500';
figureDir = fullfile(dataDir,'Figures\');
% Loading the necessary files
if ~loadTriggerData(dataDir)
    fprintf(1,'Not possible to load all the necessary variables\n')
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

promptStrings = {'Viewing window (time lapse) [s]:','Response window [s]',...
    'Bin size [s]:'};
defInputs = {'0.1, 0.1', '0.002, 0.05', '0.001'};
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
[discStack, cst] = getStacks(spkLog,Conditions(chCond).Triggers,onOffStr,...
    timeLapse,fs,fs,spkSubs,continuousSignals);

[Ne, Nt, NTa] = size(discStack);
% Computing the time axis for the stack
tx = (0:Nt - 1)/fs - timeLapse(1);
%% LTP Exclusive: Control and post-induction conditions
% Boolean flags indicating which trigger belongs to which condition (delay
% flags)
% Determining if the gap between stimuli are big enough to separate the
% conditions. If the gap is bigger than 3 standard deviations, the gap is
% sufficiently big to account for another set
[biGaps, gapSubs] = sort(diff(Conditions(chCond).Triggers(:,1)), 'descend');
gapFlag = abs(zscore(biGaps)) > 3;
mainGaps = sort(gapSubs(gapFlag), 'ascend');
Ng = numel(mainGaps);
Ncond = Ng+1;
condFlags = false(NTa, Ncond);
initSub = 1;
for cg = 1:Ng
    condFlags(initSub:mainGaps(cg),cg) = true;
    initSub = mainGaps(cg) + 1;
end
condFlags(initSub:NTa, Ncond) = true;
Na = sum(condFlags,1);
% delayFlags = condFlags;


%% Counting spikes in given windows and computing the statistical significance
% Time logical flags
sponActStackIdx = tx >= spontaneousWindow(1) & tx <= spontaneousWindow(2);
respActStackIdx = tx >= responseWindow(1) & tx <= responseWindow(2);
timeFlags = [sponActStackIdx;respActStackIdx];
% Time window
delta_t = diff(responseWindow);
% Statistical tests
[Results, Counts] = statTests(discStack,condFlags,timeFlags);
% Firing rate for all clusters, for all trials

% Plotting
figs = scatterSignificance(Results, Counts, condNames, delta_t, sortedData(goods,1));
meanfr = cellfun(@(x) mean(x,2)/delta_t,Counts,'UniformOutput',false);
%% Plots
mxfr = cellfun(@(x) max(x),meanfr);
mxfr = round(mxfr*1.15, -1);
Nr = numel(Results);
figs = gobjects(Nr,1);
ax = gobjects(2,1);
aslSubX = 1;
aslSubY = 2;
axLabels = {'Spontaneous_', 'Evoked_'};
Ncond = size(condFlags,2);
Hc = false(Ne-1, Nr*2 - Nccond);
hCount = 1;
for cr = 1:Nr
    combCell = textscan(Results(cr).Combination,'%d %d\t%s');    
    cond1 = double(combCell{1}); cond2 = double(combCell{2});
    figs(cr) = figure('Color',[1,1,1],'Visible','off');
    actvty = Results(cr).Activity(1).Type;
    if contains(actvty,'condition')
        figType = 1;
    else
        figType = 2;
    end
    csp = 1;
    while csp <= figType
        actvty = Results(cr).Activity(csp).Type;
        H = Results(cr).Activity(csp).Pvalues < 0.05;
        Hc(:,hCount) = H;
        hCount = hCount + 1;
        ttle = sprintf('%s: %d vs. %d',actvty,cond1,cond2);
        ax(csp) = subplot(1,figType,csp,'Parent',figs(cr));
        switch actvty
            case 'Spontaneous'
                aslSubX = 1; aslSubY = 1;
            case 'Evoked'
                aslSubX = 2; aslSubY = 2;
            otherwise
                aslSubX = 1; aslSubY = 2;
        end
        xaxis = meanfr{cond1, aslSubX}; yaxis = meanfr{cond2, aslSubY};
        scatter(ax(csp),xaxis,yaxis); grid(ax(csp), 'on');
        text(ax(csp), xaxis, yaxis, sortedData(goods,1))
        grid(ax(csp), 'minor'); axis('square'); 
        axis(ax(csp), [0, mxfr(cond1,aslSubX), 0, mxfr(cond1,aslSubY)]);
        ax(csp).NextPlot = 'add';
        axMx = max(mxfr(cond1,aslSubX),mxfr(cond2,aslSubY));
        line(ax(csp), 'XData', [0, axMx],...
            'YData', [0, axMx],...
            'Color', [0.8, 0.8, 0.8], 'LineStyle', '--');
        title(ax(csp),ttle); 
        xlabel(ax(csp), [axLabels{aslSubX},num2str(cond1),' [Hz]']); 
        ylabel(ax(csp), [axLabels{aslSubY},num2str(cond2),' [Hz]']);
        scatter(ax(csp),xaxis(H), yaxis(H),15, 'DisplayName', 'Significant')
        [mdl,~,rsq] = fit_poly(xaxis, yaxis, 1);
        line(ax(csp),'XData',[0, axMx],...
            'YData', [0*mdl(1) + mdl(2), axMx*mdl(1) + mdl(2)],...
            'DisplayName', 'Trend line', 'LineStyle', ':', 'LineWidth', 0.5)
        csp = csp + 1;
        if aslSubX == 1 && aslSubY == 2
            H2 = Results(cr).Activity(csp).Pvalues < 0.05;
            scatter(ax(csp-1), xaxis(H2), yaxis(H2), '.', 'DisplayName', 'Shuffled')
        end 
    end
end
set(figs,'Visible','on')

%% Waveforms
wfclIdx = any(Hc, 2);
clWf = getClusterWaveform(sortedData(goods(wfclIdx),1), dataDir);

%% Plotting the scatter data points
% % Auxiliary variable for labelling the data points in the scatter plots.
% clsLbls = num2str((1:Ncl)');
% llx = round(llx*1.1,-1);
% % Rate plots
% controlFig = figure('Color',[1,1,1],'Name','Control condition','Visible','off');
% postFig = figure('Color',[1,1,1],'Name','Post-induction','Visible','off');
% sponFig = figure('Color',[1,1,1],'Name','Spontaneous','Visible','off');
% respFig = figure('Color',[1,1,1],'Name','Evoked','Visible','off');
% controlAx = axes('Parent',controlFig);
% postAx = axes('Parent',postFig);
% sponAx = axes('Parent',sponFig);
% respAx = axes('Parent',respFig);
% 
% % Control plot (spontaneous vs evoked)
% scatter(controlAx, contRateSpontan, contRateResponse, 'DisplayName',...
%     'Rate per cluster'); grid(controlAx, 'on');
% xlabel(controlAx, 'Spontaneous rate [Hz]');
% ylabel(controlAx, 'Evoked rate [Hz]');
% title(controlAx, 'Control condition: spontaneous vs evoked')
% hold(controlAx,'on'); scatter(controlAx, contRateSpontan(HContCond == 1),...
%     contRateResponse(HContCond == 1), 'Marker','square',...
%     'DisplayName','Sync')
% scatter(controlAx, contRateSpontan(HShufContCond == 1),...
%     contRateResponse(HShufContCond == 1), 'Marker','+',...
%     'DisplayName','Shuf')
% legend(controlAx,'show'); axis(controlAx, 'square');  
% line(controlAx, 'XData', [0, llx], 'YData', [0, llx], 'Color', [0.8, 0.8, 0.8],...
%     'LineStyle', '--');
% text(controlAx, contRateSpontan, contRateResponse, clsLbls);
% controlFig.Visible = 'on';
% % Post-induction (spontaneous vs evoked)
% scatter(postAx, postRateSpontan, postRateResponse, 'DisplayName',...
%     'Rate per cluster'); grid(postAx, 'on');
% hold(postAx,'on'); scatter(postAx, postRateSpontan(HPostCond == 1),...
%     postRateResponse(HPostCond == 1), 'Marker','square',...
%     'DisplayName','Sync')
% scatter(postAx, postRateSpontan(HShufPostCond == 1),...
%     postRateResponse(HShufPostCond == 1), 'Marker','+',...
%     'DisplayName','Shuf')
% xlabel(postAx, 'Spontaneous rate [Hz]');
% ylabel(postAx, 'Evoked rate [Hz]');
% title(postAx, 'Post-induction condition: spontaneous vs evoked')
% legend(postAx,'show');  axis(postAx, 'square');  
% line(postAx, 'XData', [0, llx], 'YData', [0, llx], 'Color', [0.8, 0.8, 0.8],...
%     'LineStyle', '--');
% text(postAx, postRateSpontan, postRateResponse, clsLbls);
% postFig.Visible = 'on';
% % Spontaneous activity (Control vs post-induction)
% scatter(sponAx, contRateSpontan, postRateSpontan); grid(sponAx, 'on');
% xlabel(sponAx, 'Control spontaneous rate [Hz]');
% ylabel(sponAx, 'Post-induction spontaneous rate [Hz]');
% title(sponAx, 'Spontaneous activity: control vs post-induction')
% hold(sponAx,'on'); scatter(sponAx, contRateSpontan(HSpon == 1),...
%     postRateSpontan(HSpon == 1),'Marker','+','DisplayName','H=1');
% legend(sponAx, 'show'); axis(sponAx, 'square');  
% line(sponAx, 'XData', [0, llx], 'YData', [0, llx], 'Color', [0.8, 0.8, 0.8],...
%     'LineStyle', '--');
% text(sponAx, contRateSpontan, postRateSpontan, clsLbls);
% sponFig.Visible = 'on';
% % Evoked activity (Control vs post-induction)
% scatter(respAx, contRateResponse, postRateResponse); grid(respAx, 'on');
% xlabel(respAx, 'Control evoked rate [Hz]');
% ylabel(respAx, 'Post-induction evoked rate [Hz]');
% title(respAx, 'Evoked activity: control vs post-induction')
% hold(respAx,'on'); scatter(respAx, contRateResponse(HResp == 1),...
%     postRateResponse(HResp == 1),'Marker','+','DisplayName','H=1');
% legend(respAx, 'show'); axis(respAx, 'square');  
% line(respAx, 'XData', [0, llx], 'YData', [0, llx], 'Color', [0.8, 0.8, 0.8],...
%     'LineStyle', '--');
% text(respAx, contRateResponse, postRateResponse, clsLbls);
% respFig.Visible = 'on';
% 
% % Probability plots
% controlFigp = figure('Color',[1,1,1],'Name','Control condition','Visible','off');
% postFigp = figure('Color',[1,1,1],'Name','Post-induction','Visible','off');
% sponFigp = figure('Color',[1,1,1],'Name','Spontaneous','Visible','off');
% respFigp = figure('Color',[1,1,1],'Name','Evoked','Visible','off');
% controlAxp = axes('Parent',controlFigp);
% postAxp = axes('Parent',postFigp);
% sponAxp = axes('Parent',sponFigp);
% respAxp = axes('Parent',respFigp);
% 
% 
% % Control plot (spontaneous vs evoked)
% scatter(controlAxp, contProbSpontan, contProbResponse, 'DisplayName',...
%     'p(s) per cluster'); grid(controlAxp, 'on');
% xlabel(controlAxp, 'p(Spontaneous)');
% ylabel(controlAxp, 'p(Evoked)');
% title(controlAxp, 'Control condition_{p(s)}: spontaneous vs evoked')
% %{ 
% % STATISTICAL SIGNIFICANCE HIGHLIGHTING
% % hold(controlAxp,'on'); scatter(controlAxp, contProbSpontan(HContCond == 1),...
% %     contProbResponse(HContCond == 1), 'Marker','square',...
% %     'DisplayName','Sync')
% % scatter(controlAxp, contProbSpontan(HShufContCond == 1),...
% %     contProbResponse(HShufContCond == 1), 'Marker','+',...
% %     'DisplayName','Shuf')
% % legend(controlAxp,'show')
% %}
% text(controlAxp, contProbSpontan, contProbResponse, clsLbls);
% axis(controlAxp, 'square');  
% line(controlAxp, 'XData', [0, 1], 'YData', [0, 1], 'Color', [0.8, 0.8, 0.8],...
%     'LineStyle', '--');
% controlFigp.Visible = 'on';
% % Post-induction (spontaneous vs evoked)
% scatter(postAxp, postProbSpontan, postProbResponse, 'DisplayName',...
%     'p(s) per cluster'); grid(postAxp, 'on');
% %{
% STATISTICAL SIGNIFICANCE HIGHLIGHTING
% hold(postAxp,'on'); scatter(postAxp, postProbSpontan(HPostCond == 1),...
%     postProbResponse(HPostCond == 1), 'Marker','square',...
%     'DisplayName','Sync')
% scatter(postAxp, postProbSpontan(HShufPostCond == 1),...
%     postProbResponse(HShufPostCond == 1), 'Marker','+',...
%     'DisplayName','Shuf')
% legend(postAxp,'show')
% %}
% xlabel(postAxp, 'p(Spontaneous)');
% ylabel(postAxp, 'p(Evoked)');
% title(postAxp, 'Post-induction condition_{p(s)}: spontaneous vs evoked')
% text(postAxp, postProbSpontan, postProbResponse, clsLbls);
% axis(postAxp, 'square');  
% line(postAxp, 'XData', [0, 1], 'YData', [0, 1], 'Color', [0.8, 0.8, 0.8],...
%     'LineStyle', '--');
% postFigp.Visible = 'on';
% % Spontaneous activity (Control vs post-induction)
% scatter(sponAxp, contProbSpontan, postProbSpontan); grid(sponAxp, 'on');
% xlabel(sponAxp, 'p(Spontaneous | Control)');
% ylabel(sponAxp, 'p(Spontaneous | Post-induction)');
% title(sponAxp, 'Spontaneous activity_{p(s)}: control vs post-induction')
% %{
% STATISTICAL SIGNIFICANCE HIGHLIGHTING
% hold(sponAxp,'on'); scatter(sponAxp, contProbSpontan(HSpon == 1),...
%     postProbSpontan(HSpon == 1),'Marker','+','DisplayName','H=1');
% %}    
% text(sponAxp, contProbSpontan, postProbSpontan, clsLbls);
% axis(sponAxp, 'square');  
% line(sponAxp, 'XData', [0, 1], 'YData', [0, 1], 'Color', [0.8, 0.8, 0.8],...
%     'LineStyle', '--');
% sponFigp.Visible = 'on';
% % Evoked activity (Control vs post-induction)
% scatter(respAxp, contProbResponse, postProbResponse); grid(respAxp, 'on');
% xlabel(respAxp, 'p(Evoked | Control)');
% ylabel(respAxp, 'p(Evoked | Post-induction)');
% title(respAxp, 'Evoked activity_{p(s)}: control vs post-induction')
% %{
% STATISTICAL SIGNIFICANCE HIGHLIGHTING
% hold(respAxp,'on'); scatter(respAxp, contProbResponse(HResp == 1),...
%     postProbResponse(HResp == 1),'Marker','+','DisplayName','H=1');
% %}
% text(respAxp, contProbResponse, postProbResponse, clsLbls);
% axis(respAxp, 'square');  
% line(respAxp, 'XData', [0, 1], 'YData', [0, 1], 'Color', [0.8, 0.8, 0.8],...
%     'LineStyle', '--');
% respFigp.Visible = 'on';
% %% Saving the results (Figures and configuration)
% pause(20);
% saveFigs =...
%     questdlg('Do you wish to save the figures?','Save','Yes','No','Yes');
% saveFlag = false;
% if strcmpi(saveFigs, 'yes')
%     saveFlag = true;
%     if ~exist(figureDir,'dir')
%         if ~mkdir(figureDir)
%             fprintf(1,'There was a problem creating the ''Figures'' folder\n');
%             fprintf(1,'Please verify the location\n')
%             saveFlag = false;
%         end
%     end
% else
%     closeFigs = questdlg('Would you like to close all figures?',...
%         'Close all', 'Yes', 'No', 'No');
%     if strcmpi(closeFigs,'yes')
%         close all
%     end
% end
% 
% saveConfig = questdlg('What about the configuration for the analysis?',...
%     'Save Configuration', 'Yes', 'No', 'Yes');
% if strcmpi(saveConfig,'yes')
%     configStruct = struct('TimeLapse',timeLapse,...
%         'ResponseWindow',responseWindow,...
%         'BinSize',binSz, 'Condition', Conditions(chCond).name);
%     configFile = [expSubfix,'_ConfigStruct.mat'];
%     if ~exist(configFile,'file')
%         save([expSubfix,'_ConfigStruct.mat'],'configStruct')
%     else
%         fprintf(1,'The configuration file exists! No configuration saved\n')
%     end
% end
% 
% if saveFlag
%     figArray = [controlFig,postFig,sponFig,respFig,...
%         controlFigp,postFigp,sponFigp,respFigp];
%     figNameCell = {'Control condition', 'Post-induction condition',...
%         'Spontaneous activity', 'Evoked activity'};
%     figNameCell = [figNameCell, ...
%         cellfun(@(x) [x,' p(s)'],figNameCell,'UniformOutput',false)];
%     arrayfun(@configureFigureToPDF, figArray);
%     for cf = 1:numel(figArray)
%         emfFile = fullfile(figureDir,[figNameCell{cf},'.emf']);
%         pdfFile = fullfile(figureDir,[figNameCell{cf},'.pdf']);
%         figFile = fullfile(figureDir,[figNameCell{cf},'.fig']);
%         fileExist = [exist(emfFile,'file'),...
%             exist(pdfFile,'file'), exist(figFile,'file')];
%         if any(fileExist)
%             ovrAns =...
%                 questdlg('Some figures exist already. Do you wish to overwrite?',...
%                 'Overwrite','Yes', 'No', 'No');
%             if strcmp(ovrAns,'No')                            
%                 fprintf(1,'No figures saved!\n')
%                 break;
%             end               
%         end
%         print(figArray(cf), emfFile, '-dmeta')
%         print(figArray(cf), pdfFile, '-dpdf', '-fillpage')
%         savefig(figArray(cf), figFile)
%     end
% end
% 
% %% Verification PSTH
% consideredConditions = [1,3];
% Nccond = 2;
% goodsIdx = ~badsIdx';
% for ccond = 1:Nccond
%     [PSTH, trig, sweeps] = getPSTH(... discStack([true;whiskerResponsiveUnitsIdx;true],:,:),timeLapse,...
%         discStack,timeLapse,...
%         ~condFlags(:,ccond),binSz,fs);
%     fig = plotClusterReactivity(PSTH,trig,sweeps,timeLapse,binSz,...
%         [{Conditions(consideredConditions(ccond)).name};... sortedData(goods(whiskerResponsiveUnitsIdx),1);{'Laser'}],...
%         sortedData(goods,1)],...
%         strrep(expName,'_','\_'));
%     configureFigureToPDF(fig);
%     print(fig,fullfile(figureDir,sprintf('%s %s.pdf',...
%         expName, Conditions(consideredConditions(ccond)).name)),...
%         '-dpdf','-fillpage')
%     print(fig,fullfile(figureDir,sprintf('%s %s.emf',...
%         expName, Conditions(consideredConditions(ccond)).name)),...
%         '-dmeta')
% end
