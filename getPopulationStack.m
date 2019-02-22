function [dPopStruct, cPopStruct, conditionStruct, configStruct] =...
    getPopulationStack(EphysPath)
%GETPOPULATIONSTACK returns the population stack from all the detected
%experiments.
% Emilio Isaías-Camacho @GrohLab 2019
dPopExp = struct('Stack',cell(1,1),'SignalIDs',cell(1,1));
%% New configuration or load file
answ = questdlg(...
    'Do you want to start a new analysis or load a previous configuration',...
    'Analysis configuration','New Analysis','Load previous','Cancel',...
    'New Analysis');
switch answ
    case 'New Analysis'
        fprintf('Please attend to the prompting messages!\n')
        [configStruct, expFileNames] = createConfigStruct(EphysPath);
        if isempty(configStruct)
            fprintf('Woah! Something went wrong while creating the ')
            fprintf('configuration structure!\nPlease, try again later.\n')
            return
        end 
        fprintf('Good job! Moving on...\n')
    case 'Load previous'
        [fName, fDir] = uigetfile('*.gcf','Select an analysis file',EphysPath);
        fprintf('Loading the configuration file... \n')
        configStruct = loadPopConfigFile(fullfile(fDir,fName));
        if isempty(configStruct)
            fprintf('Woah! Something went wrong while loading the ')
            fprintf('configuration structure!\nPlease, try again later.\n')
            return
        end
        fprintf('done!\n')
        anMatDir = fullfile(EphysPath,'EphysData','AnalysisMatFiles');
        dbFileName = fullfile(EphysPath,'ephys_database.mat');
        load(dbFileName,'RecDB')
        Nx = size(RecDB,1);
        fFlag = false(Nx,1);
        for cct = 1:numel(configStruct.CellType)
            fFlag(RecDB.PhysioNucleus == configStruct.CellType{cct}) = true;
        end
        expNames = RecDB.Properties.RowNames(fFlag);
        expFileNames = cellfun(@fullfile,...
            repmat({anMatDir},sum(fFlag),1),...
            expNames,...
            cellfun(@strcat,expNames,...
            repmat({'analysis.mat'},sum(fFlag),1), 'UniformOutput',false),...
            'UniformOutput',false);
    otherwise
        fprintf('See you next time!\n')
        return
end
%% Get experiment stacks
fprintf('Building the stacks\n')
Nexp = numel(expFileNames);
dPopExp = repmat(dPopExp,Nexp,1);
cPopExp = dPopExp;
numDSig = zeros(Nexp,1,'single');
numCSig = numDSig;
Naps = numDSig;
fsExp = zeros(Nexp,2,'single');
% Cut experiment
for cexp = 1:Nexp
    fprintf('---------- Experiment Stack -----------\n')
    fprintf('%s\n',expFileNames{cexp})
    [dPopExp(cexp), cPopExp(cexp), fsExp(cexp,:)] =...
        signal_creTriggerase(configStruct,expFileNames{cexp});
    numDSig(cexp) = numel(dPopExp(cexp).SignalIDs);
    numCSig(cexp) = numel(cPopExp(cexp).SignalIDs);
    Naps(cexp) = size(dPopExp(cexp).Stack,3);
end
fs = mean(fsExp(:,1));
fsC = mean(fsExp(:,2));
Naps = [0;Naps];
Nt = size(dPopExp(1).Stack,2);
Ns = unique(numDSig);
Na = sum(Naps);
NtC = size(cPopExp(1).Stack,2);
NsC = unique(numCSig);
if numel(Ns) == 1
    % Allocate stack memory
    dPopStack = false(Ns,Nt,Na);
    cPopStack = zeros(NsC,NtC,Na,'single');
    % Compile stacks into one population stack
    for cexp = 1:Nexp
        stcSubs = sum(Naps(1:cexp)) + 1:sum(Naps(1:cexp+1));
        % If it is not the first assignment
        if ~(cexp == 1)
            % If the sinal ID do not match with the previous
            if ~(sum(strcmpi(dPopExp(cexp-1).SignalIDs,...
                    dPopExp(cexp).SignalIDs)) == Ns)
                % Re-order it
                fprintf('The order of the current signals disagrees with')
                fprintf(' the previous experiment.\nTrying to re-order it')
                fprintf('...\n')
                % If the re-ordering process did not work
                if false
                    fprintf('Bad news: Couldn''t re-order it.\n') %#ok<UNRCH>
                    fprintf('Pretending to write information in log file\n')
                    % Write in log
                    continue 
                end
            end
        end
        dPopStack(:,:,stcSubs) = dPopExp(cexp).Stack;
        cPopStack(:,:,stcSubs) = cPopExp(cexp).Stack;
    end
else
    fprintf('Haven''t implemented this possibility yet...\n')
end
fprintf('Finished building the stacks\n')
%% Exclusion of the undesired variables
fprintf('Refining the stacks: excluding and conditioning...\n')
xIdx = false(numel(dPopExp(1).SignalIDs),1);
for cev = 1:numel(configStruct.Exclude)
    xIdx = xIdx | strcmpi(dPopExp(1).SignalIDs,configStruct.Exclude{cev});
end
aux1 = sum(dPopStack(xIdx,:,:),2); % Sum up across time within the excluding variables
aux1 = squeeze(aux1); % Reshaping the stack into 2 dimensions (excludeVars x trials)
aux1 = sum(aux1,1); % Summing up across excluding variables
excludeIdx = aux1 > 0; % Exclude all trial which has a single occurance of any excluding variables
%% Conditioning the experimental stacks
% Number of conditioning variables: Ncv
Ncv = numel(configStruct.ConditionWindow);
cvIdx = false(Ncv,Na);
if ~strcmpi(configStruct.ConditionWindow(1).Name,'none')
    cwCell = squeeze(struct2cell(configStruct.ConditionWindow(:)));
    cwSubs = cell2mat(cwCell(2,:)');
    for ccv = 1:Ncv
        conWin = cwSubs(ccv,:);
        conSubs = (conWin + configStruct.ViewWindow(1))*fs + 1;
        ccvIdx = strcmpi(dPopExp(1).SignalIDs,cwCell{1,ccv});
        if conWin(1) < conWin(2)
            % The most expected result of the conditioning windows.
            cvIdx(ccv,:) = squeeze(sum(dPopStack(...
                ccvIdx,conSubs(1):conSubs(2),:),2)) > 0;
        elseif conWin(2) < conWin(1)
            % Interesting play of the times.
            timeIdx = false(1,Nt);
            timeIdx([1:conSubs(2),conSubs(1):Nt]) = true;
            cvIdx(ccv,:) = squeeze(sum(dPopStack(...
                ccvIdx,timeIdx,:),2)) > 0;
            % There's no conditioning variable.
        else
            % This is maybe the most useless selection of the conditioning
            % windows. t_1 = t_2
            fprintf('t_1 is equal to t_2. Focusing on one sample point in')
            fprintf(' the whole trial seems absurd and inaccurate.\n')
        end
    end
end
%% Wrapping up the output
dPopStruct = struct('Stack',dPopStack,'SignalIDs',{dPopExp(1).SignalIDs});
cPopStruct = struct('Stack',cPopStack,'SignalIDs',{cPopExp(1).SignalIDs});
conditionStruct = struct('ExcludeFlags',excludeIdx,'CVDFlags',cvIdx);
configStruct.SamplingFrequencies = mean(fsExp,1);
fprintf('The results are ready.\n')
end