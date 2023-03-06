function [nwbObj] = getSessionInfo4NWB(sessionPath)

p = inputParser;
addRequired(p, 'sessionPath', @(x) exist(x, 'dir'))
parse(p, sessionPath)
sessionPath = p.Results.sessionPath;

pathHere = @(x) fullfile(sessionPath, x);
look4this = @(x) dir(pathHere(x));

pathParts = strsplit(sessionPath, filesep);
batchFlag = contains(pathParts, "Batch");
pathParts = pathParts(find(batchFlag):end);
ephFlag = regexp(pathParts, "ephys\w*", "match");

rFiles = look4this("Recording*.bin");
dt = arrayfun(@(x) getDates(x, "Recording",".bin"), rFiles); 
dt = sort(dt); dt.Format = "default";

nwbObj = NwbFile();
nwbObj.session_start_time = dt(1);
nwbObj.identifier = '';
nwbObj.general_lab = 'Prof. Alexander Groh';
nwbObj.general_institution = ['Medical Biophysics,'...
    ' Institute for Physiology and Pathophysiology, Heidelberg University'...
    ', Germany.'];
nwbObj.session_description = sprintf("Awake head-fixed animal `%s` " + ...
    "on the roller, session `%s`", animal, session);
end