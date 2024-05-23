function [behStruct] = getWhiskANoseFromTable(BodyParts_angle)
%GETWHISKANOSEFROMTABLE computes the mean whisker position for both sides
%from 4 different followed whiskers.
% Emilio Isaías-Camacho @GrohLab 2022 v1 2024 v2

rw_bodyparts = cellfun(@(c) ~isempty(c), regexp( ...
    BodyParts_angle.Properties.VariableNames, 'rw\d' ) );
lw_bodyparts = cellfun(@(c) ~isempty(c), regexp( ...
    BodyParts_angle.Properties.VariableNames, 'lw\d' ) );

% Left whiskers - Ipsilateral whiskers - Stimulated whiskers
lw_mu = mean( BodyParts_angle{:,lw_bodyparts}, 2 );
lw_arc = abs( diff( BodyParts_angle{:, lw_bodyparts}(:, [1,end]), 1, 2 ) );

% Right whiskers - Contralateral whiskers - Non-stimulated whisker
rw_mu = mean( BodyParts_angle{:,rw_bodyparts}, 2 );
rw_arc = abs( diff( BodyParts_angle{:, rw_bodyparts}(:, [1,end]), 1, 2 ) );

% Whiskers-reach arc
warc = sum( 180 - [lw_mu, rw_mu], 2 );

% Whisker symmetry: positive values will be loaded towards the left side
% (puff side) and negative toward the right side (non-stimulated side).
% Zero is when the whiskers are symmetric.
wsymmetry = diff( cosd( [rw_mu, lw_mu] ), 1, 2 );

% Nose signal
nose = BodyParts_angle{:,"nose"};

% DLC Trigger cut
sNames = ["Stim-whiskers mean", "Stim-whisker fan arc", ...
    "Nonstim-Whiskers mean", "Nonstim-whisker fan arc", ...
    "Interwhisker arc", "Symmetry", "Nose"];

behSignals = cat(2, lw_mu, lw_arc, rw_mu, rw_arc, warc, wsymmetry, nose);
behStruct = struct('Names', sNames, 'Signals', behSignals);
end