function [PSTH, trig, sweeps, timeAxis] =...
    getPSTH(discreteStack, timeLapse, kIdx, binSz, fs)
%{ 
GETPSTH returns a peri-stimulus triggered histogram given a triggered
aligned stack in the form MxNxT, where M is the number of channels to
align, N is the number of samples (time) and T is the number of aligning
triggers found in the channel. 

[PSTH, triggerSummary, sweeps, timeaxis] =
    getPSTH(discreteStack, kickOutIdx, binSize (in seconds), samplingFrequency) 
    
INPUTS:
    discreteStack - AxNxM array which contains the information about the
    discrete events occurring at the given trials
    kickOutIdx - is a 1xA logic array which serves as flag to exclude the
    selected trials.
    binSize - is an integer to bin the discrete events in a smaller time
    window
    samplingFrequency - is the signal acquisition sampling frequency.

OUTPUTS: PSTH - is NxM array. N signals, M time points. tiggerSummary - is
a 1xM array with M time points. sweeps - is an integer containing the trial
number. timeaxis - is a 1xM array.

The kicking out process needs to be implemented in a dynamical form. The
events after the (first) spike channel are a guide to kick row of the
spikes out. For now, we will kick out all of the rows which contain a true
in the observed time.
%}
if binSz > 1
    disp('Assuming bin size given in milliseconds.')
    binSz = binSz * 1e-3;
end

% Computing the size of the PSTH stack
[Ne, Nt, Na] = size(discreteStack);
% Getting the raw PSTHs and kicking out the undesired trials
auxCounts = sum(discreteStack(:,:,~kIdx),3);
% Binning process.
binEls = ceil(binSz * fs);
PSTH = zeros(Ne-1,ceil(Nt/binEls));
% trig = zeros(1,ceil(Nt/binEls));
cb = 0;
sweeps = Na - sum(kIdx);
trig = auxCounts(1,:)/sum(~kIdx);
% Binned time axis
timeAxis = seconds(linspace(-timeLapse(1),timeLapse(2),ceil(Nt/binEls)));

for ce = 2:Ne
    while cb < Nt/binEls - 1
        PSTH(ce-1,cb+1) = sum(auxCounts(ce,cb*binEls+1:(cb+1)*binEls));
        % Question: to bin or not to bin the trigger signal?
        % trig(cb+1) = sum(auxCounts(1,cb*binEls+1:(cb+1)*binEls));
        cb = cb + 1;
    end
    cb = 0;
end
end