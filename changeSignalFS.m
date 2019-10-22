function [newSig, ntx] = changeSignalFS(sig, oldFS, newFS)
Ns = length(sig);
tx = (0:Ns-1) * 1/oldFS;
ntx = 0:1/newFS:tx(Ns);
sigTS = timeseries(sig,tx);
sigTS = sigTS.resample(ntx);
newSig = iir50NotchFilter(double(sigTS.Data), newFS);
end