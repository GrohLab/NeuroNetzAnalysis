function rollerposition = readRollerPositionsFile(filepath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fID = fopen(filepath, 'r');
rollerline = textscan(fID, '%d,2021-%*d-%*dT%d:%d:%f+02:00');
fclose(fID);
if all(cellfun(@(x) ~isempty(x), rollerline))
    % Roller file corresponding to the Bonsai timestamps
    rollerx = int16(rollerline{1}); rollert = string(rollerline{2})+":"+...
        rollerline{3}+":"+rollerline{4}; rollert = duration(rollert);
    rollert.Format = 'hh:mm:ss.SSSSSSS';
    rollerposition = table(rollerx, rollert, ...
        'VariableNames', {'RollerX','RollerT'});
else
    % Roller file corresponding to Arduino's micros() function.
    rollerposition = importdata(filepath);
    rollerposition(:,1) = unwrap(rollerposition(:,1), 2^15);
    rollerposition(:,2) = unwrap(rollerposition(:,2), 2^31);
end
end