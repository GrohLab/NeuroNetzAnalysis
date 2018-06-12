function [PSTH, trig, sweeps] = getPSTH(PSTHstack, kIdx, binSz, fs)
% GETPSTH returns a peri-stimulus triggered histogram given a triggered
% aligned stack in the form MxNxT, where M is the number of channels to
% align, N is the number of samples (time) and T is the number of aligning
% triggers found in the channel. 
% The kicking out process needs to be implemented in a dynamical form. The
% events after the (first) spike channel are a guide to kick row of the
% spikes out. For now, we will kick out all of the rows which contain a
% true in the observed time. 

% Computing the size of the PSTH stack
[~, Nt, Na] = size(PSTHstack);
auxCounts = sum(PSTHstack(:,:,~kIdx),3);
% Binning process.
binEls = round(binSz * fs);
PSTH = zeros(1,ceil(Nt/binEls));
% trig = zeros(1,ceil(Nt/binEls));
cb = 0;
sweeps = Na - sum(kIdx);
trig = auxCounts(1,:)/sum(~kIdx);
while cb < Nt/binEls - 1
    PSTH(cb+1) = sum(auxCounts(2,cb*binEls+1:(cb+1)*binEls));
%     trig(cb+1) = sum(auxCounts(1,cb*binEls+1:(cb+1)*binEls));
    cb = cb + 1;
end
end