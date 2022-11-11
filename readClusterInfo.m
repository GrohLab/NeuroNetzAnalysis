function clInfo = readClusterInfo(fileName)
clInfo = [];
if ~exist(fileName, 'file')
    fprintf(1, 'The given file doesn''t exist. Please verify it\n')
    return
end
if ~copyfile(fileName,strrep(fileName, '.tsv', '.csv'))
    fprintf(1,'There was an issue with the file.\n')
    fprintf(1,'Please, verify your permissions in the folder.\n')
    return
end
tempFileName = strrep(fileName, '.tsv', '.csv');
clInfo = readtable(tempFileName,...
    'ReadVariableNames', true, 'ReadRowNames', true, 'Delimiter', '\t');
clInfo.Properties.DimensionNames = {'Clusters', 'Measures'};
delete(tempFileName)
end