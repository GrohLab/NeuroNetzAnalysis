function [evals] = getSpectralEntropy(signal, dim)
%GETSPECTRALENTROPY calculates the entropy for the given spectrum
%   Detailed explanation goes here
if ~exist('dim','var') || isempty(dim)
    dim = 1;
end
[rs, cs] = size(signal);
Nwf = (rs*(dim==2)) + (cs*(dim==1));
evals = zeros(Nwf,1);
ft = fftshift(fft(signal,[],dim),dim);

for cs = 1:Nwf
    if dim == 1
        p = ft(:,cs);
    else
        p = ft(cs,:);
    end
    evals(cs) = getEntropyFromPDF(abs(p));
end
end

