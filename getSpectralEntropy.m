function [evals] = getSpectralEntropy(signal, dim)
%GETSPECTRALENTROPY calculates the entropy for the given spectrum
%   Detailed explanation goes here
if ~exist('dim','var') || isempty(dim)
    dim = 1;
end
mg2db = @(x) 20*log10(abs(x));
[rs, nc] = size(signal);
Nwf = [rs, nc]*([dim-1;dim]==1);
Ns = [rs, nc]*([dim;dim-1]==1);
wx = log10(ceil(Ns/2):Ns)';
posIdx = [false(floor(Ns/2),1);true(ceil(Ns/2),1)];
evals = zeros(Nwf, 1);
ft = mg2db(fftshift(fft(signal,[],dim),dim));
gof = evals;
mdl = zeros(2, Nwf);
for cs = 1:Nwf
    if dim == 1
        p = ft(:,cs);
    else
        p = ft(cs,:);
    end
    evals(cs) = getEntropyFromPDF(p);
    mdl(:,cs) = fit_poly(wx, p(posIdx), 1);
    gof(cs) = goodnessFit([wx, p(posIdx)], mdl(:,cs));
end
end