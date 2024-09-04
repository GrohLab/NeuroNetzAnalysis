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
my_xor = @(x) xor( x(:,1), x(:,2) );
if binSz > 1
    disp('Assuming bin size given in milliseconds and not in seconds.')
    binSz = binSz * 1e-3;
end

% Computing the size of the PSTH stack
[Ne, Nt, Na] = size(discreteStack);
% Getting the raw PSTHs and kicking out the undesired trials
%auxCounts = sum(discreteStack(:,:,~kIdx),3);
Nb = ceil(diff(timeLapse)/binSz);
sweeps = Na - sum(kIdx);
trig = sum(discreteStack(1,:,~kIdx),3)/sum(~kIdx);

txfs = (0:Nt-1)/fs + timeLapse(1);
condStack = discreteStack(2:end,:,~kIdx);
PSTH = zeros( Ne-1, Nb );
bin_centres = (1:size(PSTH,2))'.^[1,0] * [1; -0.5] * binSz + timeLapse(1);
bin_edges = [bin_centres - binSz/2, bin_centres + binSz/2];
parfor cb = 1:Nb
    bidx = my_xor( (txfs < bin_edges(cb,:)')' );
    PSTH(:,cb) = sum( condStack( :, bidx, : ), [2,3] );
end
% for ce = 2:Ne
%     tmVals = arrayfun(@(x,y) repmat(x,y,1), txfs, auxCounts(ce,:),...
%         'UniformOutput',0); tmVals = cat(1,tmVals{:});
%     h = histcounts(tmVals,'BinWidth',binSz,'BinLimits',timeLapse);
%     try
%         PSTH(ce-1,:) = h;
%     catch
%         PSTH(ce-1,1:numel(h)) = h;
%     end
% end
end