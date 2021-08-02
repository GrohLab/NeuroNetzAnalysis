function [stableCluster, firingMdl] = getMatrixStabilityMeasure(countMatrix, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Parse inputs
p = inputParser;

defSlope = 4e-3;
checkSlope = @(x) all([isnumeric(x), numel(x) == 1, x > 0]);

p.addRequired('countMatrix', @isnumeric);
p.addParameter('slopeRange', defSlope, checkSlope);

p.parse(countMatrix, varargin{:});

countMatrix = p.Results.countMatrix;
slopeRange = p.Results.slopeRange;

%% Function
[Ncl, Nt] = size(countMatrix);
%vmr = arrayfun(@(x) var(countMatrix(x,:), [], 2)./mean(countMatrix(x,:), 2),...
%    (1:size(countMatrix,1))');
firingMdl = arrayfun(@(x) fit_poly(1:Nt, countMatrix(x,:), 1), (1:Ncl)', 'UniformOutput', 0);
firingMdl = cat(2,firingMdl{:}); firingMdl = firingMdl';
% stableCluster = arrayfun(@(x,y) all([x >= 1, abs(y) < 5e-3]), vmr, firingMdl(:,1));
stableCluster = arrayfun(@(y) abs(y) < slopeRange, firingMdl(:,1));
end