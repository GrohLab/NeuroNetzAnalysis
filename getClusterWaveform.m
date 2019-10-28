function waveform = getClusterWaveform(clusterID)
%GETCLUSTERWAVEFORM reads the raw binary file and compiles the voltage
%traces for the given cluster. The output is a cell (or a structure)
%containing the mean waveform and its standard deviation. 
%   waveform = getClusterWaveform(clusterID)
%       INPUTS
%           - clusterID - character array, double or cell array containing
%           the clusters from which the waveforms are required.
%       OUTPUT
%           - waveform - cell array or double vector containing the mean
%           wavefrom for the given cluster(s).
% Emilio Isaias-Camacho @GrohLab 2019

% 

waveform = [];
checkNature = @(x) [ischar(x), iscell(x), isnumeric(x)];
if ~all(checkNature(clusterID))
    fprintf(1,'Unsupported input!\n')
    fprintf(1,'Be sure you input either the cluster ID as a ')
    fprintf(1,'character vector, a number or\nas a cell array')
    return
end
switch b12de(checkNature(clusterID))
    case 1
        fprintf(1,'Numeric ID detected\n')
    case 2
        fprintf(1,'Character ID detected\n')
    case 4
        fprintf(1,'Cell ID detected\n')
        cellfun(@ischar, clusterID)
end

end

