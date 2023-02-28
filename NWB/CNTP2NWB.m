function [outputArg1,outputArg2] = CNTP2NWB(channelMapPath, varargin)
%CNTP2NWB creates electrode 
%   Detailed explanation goes here
p = inputparser; p.CaseSensitive = false; p.KeepUnmatched = true;

validateString = @(x) isstring(x) | ischar(x);
validatePath = @(x) validateString(x) & exist(x,"file");
validateCoords = @(x) isnumeric(x) & numel(x)==2 & isvector(x);

addRequired(p, 'channelMapPath', validatePath)
addParameter(p, 'BrainStructure', 'unknown', validateString)
addParameter(p, 'Coordinates', [0,0], validateCoords)
addParameter(p, 'Company', 'Cambridge NeuroTech',validateString)
addParameter(p, 'Filtering', '', validateString)
addParameter(p, 'Comment', '', validateString)

parse(p, channelMapPath, varargin{:});

channelMapPath = p.Results.channelMapPath;
coords = p.Results.Coordinates;
brainStructure = p.Results.BrainStructure;
company = p.Results.Company;
filtering = p.Results.Filtering;    
comment = p.Results.Comment;

%fnOpts = {'UniformOutput', false};
varsInCMF = {'chanMap', 'chanMap0ind', 'connected', 'kcoords', 'name', ...
    'shankInd', 'xcoords', 'ycoords'};

if ~validatePath(channelMapPath)
    fprintf(1, "Bad input. Didn't create NWB object.\n")
    return
end
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
    'name', probeName, ...
    'description', comment, ...
    'manufacturer', company ...
);
%nwb.general_devices.set('array', device);

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