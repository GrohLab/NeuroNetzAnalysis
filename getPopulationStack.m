function [dPopStruct, cPopStruct] = getPopulationStack(EphysPath)
%GETPOPULATIONSTACK returns the population stack from all the detected
%experiments.
dPopStruct = struct('Stack',cell(1,1),'SignalIDs',cell(1,1));
cPopStruct = struct('Stack',cell(1,1),'SignalIDs',cell(1,1));
%% New configuration or load file
answ = questdlg(...
    'Do you want to start a new analysis or load a previous configuration',...
    'Analysis configuration','New Analysis','Load previous','Cancel',...
    'New Analysis');
switch answ
    case 'New Analysis'
        [configStruct, expFileNames] = createConfigStruct(EphysPath);
        if isempty(configStruct)
            fprintf('Woah! Something went wrong while creating the ')
            fprintf('configuration structure!\nPlease, try again later.\n')
            return
        end 
    case 'Load previous'
        [fName, fDir] = uigetfile('*.gcf','Select an analysis file',EphysPath);
        configStruct = loadPopConfigFile(fullfile(fDir,fName));
        if isempty(configStruct)
            fprintf('Woah! Something went wrong while loading the ')
            fprintf('configuration structure!\nPlease, try again later.\n')
            return
        end
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
Nexp = numel(expFileNames);
dPopStruct = repmat(dPopStruct,Nexp,1);
cPopStruct = dPopStruct;
numDSig = zeros(Nexp,1,'single');
numCSig = numDSig;
Naps = numDSig;
for cexp = 1:Nexp
    fprintf('---------- Experiment Stack -----------\n')
    fprintf('%s\n',expFileNames{cexp})
    [dPopStruct(cexp), cPopStruct(cexp)] =...
        signal_creTriggerase(configStruct,expFileNames{cexp});
    numDSig(cexp) = numel(dPopStruct(cexp).SignalIDs);
    numCSig(cexp) = numel(cPopStruct(cexp).SignalIDs);
    Naps(cexp) = size(dPopStruct(cexp).Stack,3);
end
%% Conditioning the experimental stacks
Nt = size(dPopStruct(1).Stack,2);
Ns = unique(numDSig);
Na = sum(Naps);
NtC = size(cPopStruct(1).Stack,2);
NsC = unique(numCSig);
if numel(Ns) == 1
    dPopStack = false(Ns,Nt,Na);
    cPopStack = false(NsC,NtC,Na);
    for cexp = 1:Nexp
        dPopStack(:,:,Naps(cexp)) = dPopStruct(cexp).Stack;
        cPopStack(:,:,Naps(cexp)) = cPopStruct(cexp).Stack;
    end
    
else
    fprintf('Haven''t implemented this possibility yet...\n')
end
fprintf('The results are ready.\n')
end