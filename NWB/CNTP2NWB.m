function [nwbObj, electrode_table, electrode_group, device] = ...
    CNTP2NWB(nwbObj, channelMapPath, varargin)
%CNTP2NWB creates 'electrodes' object for the NWB framework.
%   [nwbObj, electrodesObj, deviceObj] = CNTP2NWB(nwbObj, channelMapPath);
%    = CNTP2NWB(..., 'Name', 'Value');
% INPUTS:
%        nwbObj - NwbFile object with ideally the session information.
%        channelMapPath - char or string variable containing the path to
%                         the session channel map file fed to kilosort.
%        <name-value pairs>
%        BrainStructure - string or char variable describing the target
%                         region.
%        Coordinates - 3x1 or 1x3 numeric array containing the coordinates
%                      from the micromanipulator in micrometers. The
%                      representation is as follows: (+x = posterior, +y =
%                      inferior, +z = subject’s right). Assuming animal and
%                      electrodes facing experimenter, the given
%                      coordinates are considered to be the left-most
%                      shank.
%        Company - string or char variable containing the manufacturing
%                  company of the silicon probe.
%        Filtering - string or char variable describing the filtering
%                    information used to record the session.
% OUTPUTS:
%        electrodesObj - NWB object containing the electrodes metadata for
%                        the given session.
%        deviceObj - NWB object containing the device metadata for the
%                    session.
% Emilio Isaias-Camacho @ GrohLab 2023

%% Validation of input variables
p = inputParser; p.CaseSensitive = false; p.KeepUnmatched = true;

validateString = @(x) isstring(x) | ischar(x);
validatePath = @(x) validateString(x) & exist(x,"file");
validateNWB = @(x) isa(x, 'NwbFile');
% Coordinates are three dimensional; (ML, AP, DV); DV must always be
% negative.
validateCoords = @(x) isnumeric(x) & numel(x)==3 & isvector(x) & x(3)<=0;

addRequired(p, "nwbObj", validateNWB)
addRequired(p, 'channelMapPath', validatePath)
addParameter(p, 'BrainStructure', 'unknown', validateString)
addParameter(p, 'Coordinates', [0,0,0], validateCoords)
addParameter(p, 'Company', 'Cambridge NeuroTech',validateString)
addParameter(p, 'Filtering', 'unknown', validateString)

parse(p, nwbObj, channelMapPath, varargin{:});

nwbObj = p.Results.nwbObj;
channelMapPath = p.Results.channelMapPath;
% Coordinates as medio-lateral and anterio-posterior
coords = p.Results.Coordinates;
brainStructure = p.Results.BrainStructure;
company = p.Results.Company;
filtering = p.Results.Filtering;
%% Auxiliary functions
fnOpts = {'UniformOutput', false};
varsInCMF = {'chanMap', 'chanMap0ind', 'connected', 'kcoords', 'name', ...
    'shankInd', 'xcoords', 'ycoords'};
tblVars = {'x', 'y', 'z', 'imp', 'location', 'filtering', 'group', ...
    'label', 'rel_x', 'rel_y', 'rel_z'};
cmMF = matfile(channelMapPath);
varsFlags = contains(varsInCMF, who(cmMF), "IgnoreCase", true);
% Assuming 1 shank if shankInd wasn't specified.
nshanks = 1;
if varsFlags(6)
    shID = cmMF.shankInd; shID = shID(:);
    nshanks = unique(shID);
end
if varsFlags(5); probeName = cmMF.name; end
%% Creation of NWB objects
device = types.core.Device('description', probeName, ...
    'manufacturer', company);
nwbObj.general_devices.set('array', device);

x = cmMF.xcoords; y = cmMF.ycoords;
% (+x = posterior, +y = inferior, +z = subject’s right).
shFlags = shID == nshanks(:)';
elec_tbl = cell2table(cell(0, length(tblVars)), 'VariableNames', tblVars);
for csh = 1:nshanks
    electrode_group = types.core.ElectrodeGroup(...
        'description', ['electrode group for shank' num2str(csh)],...
        'location', brainStructure,...
        'device', types.untyped.SoftLink(device));
    nwbObj.general_extracellular_ephys.set(['shank' num2str(csh)], ...
        electrode_group);
    group_object_view = types.untyped.ObjectView(electrode_group);
    for ce = find(shFlags(:,csh))'
        electrode_label = ['shank' num2str(csh) 'elec' num2str(ce)];
        elec_tbl = [elec_tbl; ...
            {coords(1), coords(2), coords(3), NaN, brainStructure, ...
            filtering, group_object_view, electrode_label, 0, y(ce),...
            x(ce)}]; %#ok<AGROW> 
    end
end
%{
electrode_group = arrayfun(@(sh) types.core.ElectrodeGroup( ...
    'description', ['electrode group for shank' num2str(sh)], ...
    'location', brainStructure, ...
    'device', types.untyped.SoftLink(device)), nshanks);
for ce = 1:numel(nshanks)
    nwbObj.general_extracellular_ephys.set("shank"+nshanks(ce), ...
        electrode_group(ce));
end

group_object_view = arrayfun(@(e) types.untyped.ObjectView(e), ...
    electrode_group); 
shFlags = shID == nshanks(:)';
group_object_views = arrayfun(@(sh) ...
    group_object_view(shFlags(shFlags(:,sh),sh)*sh), nshanks, fnOpts{:});
group_object_views = cat(1, group_object_views{:});
% Ordering the electrode labels
eNum = arrayfun(@(s) cumsum(shFlags(shFlags(:,s),s)), nshanks, fnOpts{:});
eNum = cat(1, eNum{:}); eSubs = zeros(size(shFlags));eSubs(shFlags) = eNum;
electrode_labels = arrayfun(@(s,e) "shank"+string(s)+"elec"+string(e), ...
    shID, sum(eSubs,2));


rm = @(x) repmat(x, size(shID,1), 1);

tbl = table(eCoords(:,1), eCoords(:,2), eCoords(:,3), nan(size(eCoords,1),1), ...
    rm(brainStructure), rm(filtering), group_object_views, electrode_labels, ...
    'VariableNames',tblVars);
%}
electrode_table = util.table2nwb(elec_tbl, 'all electrodes');
nwbObj.general_extracellular_ephys_electrodes = electrode_table;
end