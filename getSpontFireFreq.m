function [firingRatePerCluster, deltaTrigTimeSum, sponSpks, sponIsi] =...
    getSpontFireFreq(spks, trigs, consTime, fs, bufTime)
%GETSPONTFIREFREQ computes the firing rate outside of the specified
%triggers leaving a time buffer around the pulse. 
%
%   [firingRatePerCluster, deltaTrigTime, sponSpks, sponIsi] = 
%                       getSpontFireFreq(spks, trigs, fs, bufTime)
%
%   INPUTS:
%       spks - spike times or subscripts for each cluster in a cell array
%       trigs - an Nx2 matrix containing all pulses that the user wants to
%               avoid counting as spontaneous
%       consTime - is a 1x2 vector containing the times enclosing the
%                  spontaneous windows
%       fs - sampling frequency of the experiment
%       bufTime - scalar indicating the time after the offset that should
%                 be omitted
%   OUTPUT:
%       firingRatePerCluster - a column vector containing the firing rate
%                              per cluster in the given order.
%       deltaTrigTimeSum - a column vector containing the time considered
%                          as spontaneous
%       sponSpks - a cell array with the spikes per cluster within the
%                  considered time in consTime and outside the given
%                  triggers in trigs.
%       sponIsi - a cell array with the inter-spike intervals for the
%                 spontaneous spikes.
% Emilio Isaias-Camacho @GrohLab 2020
fnOpts = {'UniformOutput', false};

spksFlag =...
    cellfun(@(x) round(x./fs) > consTime(1) & round(x./fs) <= consTime(2),...
    spks, fnOpts{:});
isi = cellfun(@diff, spks, fnOpts{:});
consSpks = cellfun(@(x,y) x(y), spks, spksFlag, fnOpts{:});
consIsi = cellfun(@(x,y) x(y(1:end-1)), isi, spksFlag, fnOpts{:});
compareTimes =@(x) arrayfun(@(z) any(trigs(:,1) < z &...
    trigs(:,2)+(bufTime*fs) >= z), x);
evokSpksFlags = cellfun(@(x) compareTimes(x), consSpks, fnOpts{:});
NSpkSpont = cellfun(@(x) sum(~x), evokSpksFlags);
trigs = trigs(trigs(:,2) <= round(consTime(2) * fs),:);
deltaTrigTimeSum = sum([trigs(:,1); round(consTime(2) * fs)] -...
    [round(consTime(1) * fs); trigs(:,2)])/fs;
firingRatePerCluster = NSpkSpont./deltaTrigTimeSum;
sponSpks = cellfun(@(x,y) x(~y), consSpks, evokSpksFlags, fnOpts{:});
sponIsi = cellfun(@(x,y) x(~y(1:numel(x))), consIsi, evokSpksFlags,...
    fnOpts{:});
end
