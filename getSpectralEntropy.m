function [evals, mdl, gof] = getSpectralEntropy(signal, dim, fs)
%GETSPECTRALENTROPY calculates the entropy for the given spectrum
%   Detailed explanation goes here
if ~exist('dim','var') || isempty(dim)
    dim = 1;
end
mg2db = @(x) 20*log10(abs(x));
[rs, nc] = size(signal);
Nwf = [rs, nc]*([dim-1;dim]==1);
Ns = [rs, nc]*([dim;dim-1]==1);
wx = log10((ceil(Ns/2)+(1-mod(Ns,2))):Ns)';
posIdx = [false(floor(Ns/2),1);true(ceil(Ns/2),1)];
evals = zeros(Nwf, 1);
ft = fftshift(fft(signal,[],dim),dim); dw = fs/size(fs,1);
mt = abs(ft);
pt = angle(ft);
dbft = mg2db(ft);
gof = evals;
mdl = zeros(2, Nwf);
for cs = 1:Nwf
    if dim == 1
        p = mt(:,cs);
        dp = dbft(:, cs);
    else
        p = mt(cs,:);
        dp = dbft(cs, :);
    end
    evals(cs) = abs(getEntropyFromPDF(dp(posIdx)));
    mdl(:,cs) = fit_poly(wx, dp(posIdx), 1);
    gof(cs) = goodnessFit([wx, dp(posIdx)], mdl(:,cs));
    wcp = getWaveformCriticalPoints(dp(posIdx), 1/dw);
end
end