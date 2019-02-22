function [binnedSignal] = binTime(signal,binSz,fs)
%BINTIME returns a binned signal using the bin size and the sampling
%frequency as input arguments. It reshapes the signal into packages of
%elements fitting into the specified bins and sums the values 'in-one-go'.
%   Emilio Isaías-Camacho @GrohLab 2019
binEls = ceil(binSz*fs);
N = length(signal);
binCnt = round(N/binEls);
rshTop = binCnt * binEls;
arrSig = reshape(signal(1:rshTop),binEls,binCnt);
binnedSignal = sum(arrSig)/binEls;
end