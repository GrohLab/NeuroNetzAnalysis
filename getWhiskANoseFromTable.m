function [behStruct] = getWhiskANoseFromTable(BodyParts_angle)
%GETWHISKANOSEFROMTABLE computes the mean whisker position for both sides
%from 4 different followed whiskers.
% Emilio Isaías-Camacho @GrohLab 2022
% Left whiskers - Ipsilateral whiskers
lw = mean(BodyParts_angle{:,{'lw1', 'lw2', 'lw3', 'lw4'}},2);
lw = lw - mean(lw);
% Right whiskers - Contralateral whiskers 
rw = mean(BodyParts_angle{:,{'rw1', 'rw2', 'rw3', 'rw4'}},2);
rw = rw - mean(rw);
% Nose signal
nose = BodyParts_angle{:,"nose"} - mean(BodyParts_angle{:,"nose"});

% DLC Trigger cut
sNames = ["Ipsi-Whiskers", "Contra-Whiskers", "Nose"];

behSignals = cat(2, lw, rw, nose);
behStruct = struct('Names', sNames, 'Signals', behSignals);
end