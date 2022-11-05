function [iok] = saveConfigFile(conFiNa, conFlags, sigNames)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if size(conFlags,2) == length(sigNames)
    fID = fopen(conFiNa,'w');
    % Writing the trigger
    fprintf(fID,'t:\t %s\n',sigNames{conFlags(:,1)});
    % Writing the 
else
    fprintf('Please verify the inputs to the function.\n')
    fprintf('The length of the configuration flag and the signal names ')
    fprintf('differ.\n')
    iok = -1;
end

end

