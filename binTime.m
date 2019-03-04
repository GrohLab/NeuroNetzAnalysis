function [binnedSignal, binWidth] = binTime(signal,binSz,fs)
%BINTIME returns a binned signal using the bin size and the sampling
%frequency as input arguments. It reshapes the signal into packages of
%elements fitting into the specified bins and sums the values 'in-one-go'.
%   Emilio Isaías-Camacho @GrohLab 2019
binWidth = ceil(binSz * fs);
[Nr, Nc] = size(signal);
if Nc > Nr
    signal = signal';
    N = Nc;
    Nc = Nr;
    Nr = N;
end
N = Nr;
binCnt = round(N/binWidth);
rshTop = binCnt * binWidth;
arrSig = reshape(signal(1:rshTop,:),binWidth,binCnt,Nc);
binnedSignal = reshape(sum(arrSig,1)./binWidth, binCnt, Nc);
if rshTop < N
    binnedSignal(binCnt+1,:) =...
        sum(signal(rshTop+1:N,:),1)/numel(rshTop+1:N);
end
end