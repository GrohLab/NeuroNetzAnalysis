outDir = 'Z:\SC Anatomy paper data\Roller';
% Session path
sessionPath = fullfile("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch10_ephys.e\MC\GADe60\221026_C_DV2150\ephys_E1");
% Coordinates from file name or labbook
%(+x = posterior AP, +y = inferior DV, +z = subject s right ML)
coords = [3600, 2150, 1500];
% Used probe
cmGeometry = 'E';
% Date (and time) of session start
sessionDate = datetime('221026','InputFormat','yyMMdd'); sessionDate.Format = 'yyMMdd';
% Identifier
identifier = ['WT60_' char(sessionDate)];
genotype = 'WT';

birthdate = datetime('21.4.22','InputFormat','dd.MM.yy');

mouseSource = 'IBF';
mouseSex = 'male';

selected_condition = 'Control Puff';
experimentDescription = 'Air puff stimulation on the whiskers of an awake head-fixed mouse on a roller';

nwbObj = convertSession2NWB(sessionPath, ...
    "SessionDate", sessionDate, ...
    "Coordinates", coords, ...
    "ConditionSelection", selected_condition,...
    "OutDir", outDir, ...
    "Birthdate", birthdate, ...
    "Identifier", identifier, ...
    "ExperimentDescription", experimentDescription, ...
    "MouseSex", mouseSex, ...
    "MouseSource", mouseSource, ...
    "Genotype", genotype, ...
    "chanMapGeometry", cmGeometry);