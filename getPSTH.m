function [PSTH, trig, sweeps, timeAxis] =...
    getPSTH(discreteStack, timeLapse, kIdx, binSz, fs)
%GETPSTH returns a peri-stimulus triggered histogram given a triggered
% aligned stack in the form MxNxT, where M is the number of channels to
% align, N is the number of samples (time) and T is the number of aligning
% triggers found in the channel. 
% 
% [PSTH, triggerSummary, sweeps, timeaxis] =
%     getPSTH(discreteStack, timeLapse, kIdx, binSz, fs) 
%     
% INPUTS:
%     discreteStack - NxMxA array which contains the information about the
%     discrete events occurring at the given trials.
%     timeLapse - 1x2 array containing the time before and after the
%     trigger to consider
%     kIdx - is a 1xA logic array which serves as flag to exclude the
%     selected trials.
%     binSz - is an integer to bin the discrete events in a smaller time
%     window
%     fs - is the signal acquisition sampling frequency.
% 
% OUTPUTS: 
%     PSTH - is NxM array. N signals, M time points. 
%     tiggerSummary - is a 1xM array with M time points. 
%     sweeps - is an integer containing the trial number. 
%     timeaxis - is a 1xM array.
% 
%   Emilio Isaias-Camacho @ GrohLab 2019

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