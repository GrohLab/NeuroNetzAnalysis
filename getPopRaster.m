function [RasterStruct] = getPopRaster(dst, conditStruct, configStruct)
%GETPOPRASTER creates the necessary variables to plot the time activation
%of each cell in the population and to see simultaneous activity of the
%conditioning variables.
%   Emilio Isaías-Camacho @ GrohLabs 2019
[Nv, Nts, Na] = size(dst.Stack);

RasterStruct = conditStruct;
end

