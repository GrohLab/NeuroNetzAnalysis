function clusterInfo = getClusterInfo(filename)
getContFirstCell = @(x) x{1};
clusterInfo = table;
if ~exist(filename,'file')
    fprintf(1,'The given file doesn''t exist\n')
    return
end

fID = fopen(filename,'r');
heads = getContFirstCell(textscan(fgetl(fID), '%s\t'));

while ~feof(fID)
    conts = getContFirstCell(textscan(fgetl(fID), '%s',...
        'Delimiter', '\t', 'Whitespace', ' '));
end
end
