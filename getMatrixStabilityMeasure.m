function [stableCluster] = getMatrixStabilityMeasure(countMatrix, varargin)
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

%% Getting measures of the matrix' rows
[Ncl, Nt] = size(countMatrix); fnOpts = {'UniformOutput', false};
fanoFact = @(x) var(x, [], 2, "omitnan")./mean(x, 2, "omitnan");
ffm = fanoFact(countMatrix);
firingMdl = arrayfun(@(x) fit_poly(1:Nt, countMatrix(x,:), 1), (1:Ncl)',...
    fnOpts{:}); firingMdl = cat(2,firingMdl{:}); firingMdl = firingMdl';
% stableCluster = arrayfun(@(y) abs(y) < slopeRange, firingMdl(:,1));
stableCluster = (abs(firingMdl(:,1)) < slopeRange) & (ffm <= 1);
end