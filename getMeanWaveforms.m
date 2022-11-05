function [meanWf, txwf] = getMeanWaveforms(clWaveCell, fs)
%GETMEANWAVEFORMS returns a matrix of NxM with the mean M waveforms from N
%clusters and a time vector of 1xN elements
gmw = @(x) mean(x,2);
meanWf = cell2mat(cellfun(gmw, clWaveCell(:,2), 'UniformOutput', 0)');
txwf = (0:size(meanWf,1)-1)/fs;
end