function [outputArg1,outputArg2] = CNTP2NWB(channelMapPath, varargin)
%CNTP2NWB creates 'electrodes' object for the NWB framework.
%   [electrodesObj, deviceObj] = CNTP2NWB(channelMapPath);
% INPUTS:
%        channelMapPath - char or string variable containing the path to
%                         the session channel map file fed to kilosort.
%        <name-value pairs>
%        BrainStructure - string or char variable describing the target
%                         region.
%        Coordinates - 3x1 or 1x3 numeric array containing the coordinates
%                      from the micromanipulator in micrometers. The
%                      representation is as follows: (ML, AP, DV), where ML
%                      is medio-lateral, AP is antero-posterior, and DV is
%                      dorso-ventral and must be negative.
%                      Assuming animal and electrodes facing experimenter,
%                      the given coordinates are considered to be the
%                      left-most shank.
%        Company - string or char variable containing the manufacturing
%                  company of the silicon probe.
%        Filtering - string or char variable describing the filtering
%                    information used to record the session.
%        Comment - String or char variable used to comment any observation.
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
% Coordinates are three dimensional; (ML, AP, DV); DV must always be
% negative.
validateCoords = @(x) isnumeric(x) & numel(x)==3 & isvector(x) & x(3)<=0;

addRequired(p, 'channelMapPath', validatePath)
addParameter(p, 'BrainStructure', 'unknown', validateString)
addParameter(p, 'Coordinates', [0,0,0], validateCoords)
addParameter(p, 'Company', 'Cambridge NeuroTech',validateString)
addParameter(p, 'Filtering', 'unknown', validateString)
addParameter(p, 'Comment', '', validateString)

parse(p, channelMapPath, varargin{:});

channelMapPath = p.Results.channelMapPath;
% Coordinates as medio-lateral and anterio-posterior
coords = p.Results.Coordinates;
brainStructure = p.Results.BrainStructure;
company = p.Results.Company;
filtering = p.Results.Filtering;    
comment = p.Results.Comment;
%% Auxiliary functions and 
%fnOpts = {'UniformOutput', false};
varsInCMF = {'chanMap', 'chanMap0ind', 'connected', 'kcoords', 'name', ...
    'shankInd', 'xcoords', 'ycoords'};

cmMF = matfile(channelMapPath);
varsFlags = contains(varsInCMF, who(cmMF), "IgnoreCase", true);
% Assuming 1 shank if shankInd wasn't specified.
nshanks = 1;
if varsFlags(6)
    shID = cmMF.shankInd;
    nshanks = unique(shID);
end
if varsFlags(5); probeName = cmMF.name; end

variables = {'x', 'y', 'z', 'imp', 'location', 'filtering', 'group', 'label'};
tbl = cell2table(cell(0, length(variables)), 'VariableNames', variables);
device = types.core.Device(...
    'description', probeName, ...
    'manufacturer', company...
);
%nwb.general_devices.set('array', device);
% x, y, and z represent medio-lateral, antero-posterior, and dorso-ventral
% axes respectively.
x = coords(1) + cmMF.xcoords; y = coords(2); z = coords(3) + cmMF.ycoords;
for ishank = 1:nshanks
    electrode_group = types.core.ElectrodeGroup( ...
        'description', ['electrode group for shank' num2str(ishank)], ...
        'location', brainStructure, ...
        'device', types.untyped.SoftLink(device) ...
    );
    nwb.general_extracellular_ephys.set(['shank' num2str(ishank)], electrode_group);
    group_object_view = types.untyped.ObjectView(electrode_group);
    shFlag = shID == ishank;
    nchannels_per_shank = sum(shFlag);
    for ielec = 1:nchannels_per_shank
        electrode_label = ['shank' num2str(ishank) 'elec' num2str(ielec)];
        tbl = [...
            tbl; ...
            {5.3, 1.5, 8.5, NaN, brainStructure, filtering, group_object_view, electrode_label} ...
        ];
    end
end
electrode_table = util.table2nwb(tbl, 'all electrodes');

end 