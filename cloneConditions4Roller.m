function CondRoller = cloneConditions4Roller(dataDir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Auxiliary variables and functions
dateFormStr = 'yyyy-MM-dd''T''HH_mm_ss';
outFileFormat = "RollerConditions%s.mat";
fnOpts = {"UniformOutput", false};
CondRoller = struct();

    function memFlags = findMemberships(tFlag)
        memFlags = arrayfun(@(x) ismember(Conditions(tFlag).Triggers(:,1),...
            Conditions(x).Triggers(:,1)), mixSubs, fnOpts{:});
        memFlags = [memFlags{:}];
    end


%% Validation section
% Given folder exists?
if ~exist(dataDir, "dir")
    fprintf(1, "Given directory doesn''t exist!\n")
    return
end

% Search pattern for the rest of the files
sp = fullfile(dataDir, "**");
% Does the folder have files with roller positions?
atFiles = dir(fullfile(sp, "ArduinoTriggers*.mat"));
if isempty(atFiles)
    fprintf(1, "No arduino trigger files found in:\n%s\n", dataDir)
    return
end
% Does the experiment contain intan trigger signals?
cfFiles = dir(fullfile(sp, "*analysis.mat"));
if isempty(cfFiles)
    fprintf(1, "No analysis files found in:\n%s\n", dataDir)
    return
end

% Checking for number of files
Nat = numel(atFiles);

%% Load conditions from the ephys
cfName = fullfile(cfFiles.folder, cfFiles.name);
auxVars = load(cfName, "Conditions");
if isempty(auxVars)
    fprintf(1, "Warning! Unable to load necessary variables")
    fprintf(1, " from the condition file:\n%s\n", cfFiles.name)
    return
end
Conditions = auxVars.Conditions;
condNames = arrayfun(@(x) string(x.name), Conditions);
pFlag = ismember(condNames, "PuffAll");
lFlag = ismember(condNames, "LaserAll");
mixSubs = find(~(pFlag | lFlag));
pmFlags = findMemberships(pFlag);
lmFlags = findMemberships(lFlag);

for ctf = 1:Nat
    atName = fullfile(atFiles(ctf).folder, atFiles(ctf).name);
    auxVars = load(atName);
    if isempty(auxVars)
        fprintf(1, "Warning! Unable to load necessary variables")
        fprintf(1, " from the condition file\n%s", cfFiles(ccf).name)
        return
    end
end


end

