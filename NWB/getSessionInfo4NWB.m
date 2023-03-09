function [nwbObj] = getSessionInfo4NWB(sessionPath, varargin)
%GETSESSIONINFO4NWB prepares a NwbFile object to store all information for
%the PLOS Biology SC anatomy paper.
%           nwbObj = getSessionInfo4NWB(sessionPath)
%           nwbObj = getSessionInfo4NWB(sessionPath, 'Comment', comment)
%   INPUT:
%       sessionPath - string or char variable containing the ephys data of
%       the current session to convert to NWB.
%       Name-Value pair: Comment - string or char variable adding
%       information about the session.
%   OUTPUT:
%       nwbObj - initialised NwbFile object 
% Emilio Isaías-Camacho @ GrohLab 2023

%% Validate inputs
istxt = @(x) isstring(x) | ischar(x) | iscellstr(x);

p = inputParser;
addRequired(p, 'sessionPath', @(x) exist(x, 'dir'))
addParameter(p, 'Comment', '', istxt)
addParameter(p, "Identifier", '', @(x) istxt(x) & contains(x, '_'))
addParameter(p, "SessionDate", [], @isdatetime);

parse(p, sessionPath, varargin{:})
sessionPath = p.Results.sessionPath;
comment = p.Results.Comment;
identifier = p.Results.Identifier;
sessionDate = p.Results.SessionDate;
%% Auxiliary functions and variables
pathHere = @(x) fullfile(sessionPath, x);
look4this = @(x) dir(pathHere(x));
%% Dig info out of the path
pathParts = strsplit(sessionPath, filesep);
batchFlag = contains(pathParts, "Batch");
pathParts = pathParts(find(batchFlag):end);
ephFlag = contains(pathParts, "ephys");

rFiles = look4this("Recording*.bin");
if ~isempty(rFiles) && isempty(sessionDate)
    dt = arrayfun(@(x) getDates(x, "Recording",".bin"), rFiles);
    dt = sort(dt); dt.Format = "default";
    animal = 'Mouse';
    if ephFlag(end)
        animal = regexp(pathParts, '[A-Z]+([a-z])?[0-9]+', 'match');
        animal = animal(~cellfun(@isempty, animal));
        if ~isempty(animal)
            if iscell(animal)
                animal = animal{:};
            end
            animal = animal(1); animal = char(animal);
        end
        sessionDate = datetime(dt(1),"Format", "yyMMdd");
    else
        sessionDate = regexp(pathParts, "[0-9]{6}","match");
        sessionDate(cellfun(@isempty, sessionDate)) = [];
        if iscell(sessionDate)
            sessionDate = sessionDate{:};
        end
        sessionDate = datetime(sessionDate, 'InputFormat','yyMMdd');
    end
    identifier = [char(animal), '_', char(sessionDate)];
elseif ~isempty(identifier)
    animal = extractBefore(identifier, '_'); 
    sessionDate = extractAfter(identifier, '_');
    sessionDate = datetime(sessionDate, 'InputFormat','yyMMdd');
    
else
    animal = 'Mouse'; sessionDate = datetime(); 
end
sessionDate.Format = "yyMMdd";

%% Create NWB object

nwbObj = NwbFile();
nwbObj.identifier = identifier;
nwbObj.general_lab = 'Prof. Alexander Groh';
nwbObj.general_experiment_description = 'Air puff whisker stimulation';
nwbObj.general_institution = ['Medical Biophysics,'...
    ' Institute for Physiology and Pathophysiology, Heidelberg University'...
    ', Germany.'];
nwbObj.session_description = sprintf(['Animal `%s`, session `%s`. ',...
    comment], animal, string(sessionDate));
nwbObj.session_start_time = sessionDate;
nwbObj.general_related_publications = [
    'Martín-Cortecero, J., Isaías-Camacho, E. U., Boztepe, B.,', ...
    ' Ziegler, K., Mease, R. A.,; Groh, A. (2023). ', ...
    'Monosynaptic trans-collicular pathways for sensory-motor ', ...
    'integration in the whisker system. DOI:10.1101/2022.08.30.505868'];
device = types.core.Device('description', 'Custom built laser control system');
nwbObj.general_devices.set('Laser', device);
nwbObj.general_optogenetics.set('OptogeneticStimulusSite', ...
    types.core.OptogeneticStimulusSite(...
    'excitation_lambda', 488, ...
    'location', 'Same as electrode', ...
    'device', types.untyped.SoftLink(nwbObj.general_devices.get('Laser')), ...
    'description', ['Laser, if present, delivered 200 microns away from',...
    ' electrodes']));
end