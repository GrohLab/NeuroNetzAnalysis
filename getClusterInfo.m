function clusterInfo = getClusterInfo(filename)
% GETCLUSTERINFO reads the phy generated file 'cluster_info.tsv' and
% returns a table with the variables assigned into an output table
% 'clusterInfo'. It will have size NxM, with N clusters and M variables.
%       clusterInfo = getClusterInfo(filename)
%           INPUTS
%               - filename - string indicating the relative or absolute
%               path to the 'cluster_info.tsv' file.
%           OUTPUTS
%               - clusterInfo - an NxM table with N rows corresponding to
%               each cluster and M columns corresponding to each variable
%               in the file (e.g. id, amplitude, channel, etc.)
% Emilio Isaias-Camacho @ GrohLab 2020

%% Validation of the input argument; the filename
% Checking for existence
if ~exist(filename,'file')
    fprintf(1,'The given file doesn''t exist\n')
    return
end
% Checking for validity
[~, ~, fext] = fileparts(filename);
fext = char(fext);
if ~strcmp(fext,'.tsv')
    fprintf(1,'Wrong file fed to the function!\n')
    fprintf(1,'It must be the phy-generated file ''cluster_info.tsv''\n')
    fprintf(1,'Try again with the correct file\n')
    return
end
%% Reading the file
% Open the file for reading
clusterInfo = readClusterInfo(filename);
heads = clusterInfo.Properties.VariableNames;
%% Scanning the extracted contents 
Nv = size(clusterInfo, 2);
for cv = 1:Nv
    % Proceeding slightly different for each variable
    switch heads{cv}
        case {'cluster_id','id','KSLabel','group','NeuronType','Region'}
            % Read as string if the header indicates strings
            vals = string(clusterInfo{:,cv});
        case {'firing_rate', 'fr'}
            % Read specially as the values are accompanied by the units
            if isnumeric(clusterInfo{:,cv})
                vals = clusterInfo{:,cv};
            else
                vals = cell2mat(cellfun(@(x) textscan(x,'%f spk/s'),...
                    clusterInfo{:,cv}));
            end
        otherwise
            % Read simply as a number for the rest of the variables
            vals = clusterInfo{:,cv};
    end
    % Assigning the values to the corresponding variable for all custers
    clusterInfo.(heads{cv}) = vals;
end
% Naming the rows and columns appropiately and respectively.
clusterInfo.Properties.DimensionNames = {'Clusters', 'Measures'};
clusterInfo.Properties.RowNames = clusterInfo.(heads{1});
end
