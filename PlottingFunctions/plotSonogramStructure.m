function fig = plotSonogramStructure(sonoStruct)
%PLOTSONOGRMSTRUCTURE takes the output from the SONOGRAM function and plots
%it accordingly. 
%   fig = plotSonogramStructure(sonoStruct)
%   INPUTS
%       sonoStruct - structure containing three matrices. The Fourier transform for each
%       overlapping window in the signal, the magnitude of this, and its
%       vectorial angle.
%   OUTPUTS
%       fig - figure object containing the plot of the sonogram.
% Emilio Isa?as-Camacho @GrohLab 2019
p = inputParser;
checkStruct = @(x) isstruct(x) && isfield(x,'SpectrumImage') &&...
    isfield(x, 'TimeAxis') && isfield(x,'FrequencyAxis');

addRequired(p,'sonoStruct',checkStruct);

p.KeepUnmatched = false;
parse(p,sonoStruct)

sonoStruct = p.Results.sonoStruct;

fig = figure('Color',[1,1,1]);
ax = axes('Parent',fig);
imagesc(ax, sonoStruct.TimeAxis, sonoStruct.FrequencyAxis,...
    20*log10(abs(sonoStruct.SpectrumImage))); colormap(gray);

end

