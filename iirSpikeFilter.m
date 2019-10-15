function [csig] = iirSpikeFilter(signal,sampling_frequency,cutFreq)
%IIRSPIKEFILTER returns the filtered signal according to either the user
%input in cutFreq or uses the default 600 - 9 kHz cut frequencies for the
%band pass. 
%   signal contains a broadband frequency spectra
%   sampling_frequency has a sensible name
%   cutFreq is a vector containing the low and high cut frequencies in this
%   order
[b, a] = cheby2(3,30,1/sqrt(2));
if ~exist('cutFreq','var') || isempty(cutFreq) || length(cutFreq) ~= 2 ||...
        cutFreq(1) > cutFreq(2)
    cutFreq = [600, 9e3];
end
wo = (2*mean(cutFreq))/sampling_frequency;
wc = (2*cutFreq)./sampling_frequency;
fprintf('Spike filter set at: [%.2f, %.2f]\n',cutFreq(1),cutFreq(2))
[num, den] = iirlp2bp(b, a, wo, wc);
fprintf('Filtering the signal...')
if ~isa(signal,'double')
    signal = double(signal);
end
csig = filtfilt(num, den, signal);
fprintf(' done!\n');
end

