function rollerposition = readRollerPositionsFile(filepath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fID = fopen(filepath, 'r');
rollerline = textscan(fID, '%d,2021-%d-%dT%d:%d:%f+02:00');
fclose(fID);
rollerx = int16(rollerline{1}); rollert = string(rollerline{4})+":"+...
    rollerline{5}+":"+rollerline{6}; rollert = duration(rollert);
rollert.Format = 'hh:mm:ss.SSSSSSS';
rollerposition = table(rollerx, rollert, ...
    'VariableNames', {'RollerX','RollerT'});

end

