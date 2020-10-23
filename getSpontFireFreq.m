function [firingRatePerCluster, deltaTrigTime] =...
    getSpontFireFreq(spks, trigs, consTime, fs, bufTime)
%GETSPONTFIREFREQ computes the firing rate outside of the specified
%triggers leaving a time buffer around the pulse. 
%
%   [firingRatePerCluster] = getSpontFireFreq(spks, trigs, fs, <bufCent>)
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

spksFlag =...
    cellfun(@(x) round(x./fs) > consTime(1) & round(x./fs) <= consTime(2),...
    spks, 'UniformOutput', 0);
consSpks = cellfun(@(x,y) x(y), spks, spksFlag, 'UniformOutput', 0);
compareTimes =@(x) arrayfun(@(z) any(trigs(:,1) < z & trigs(:,2)+(bufTime*fs) >= z), x);
evokSpksFlags = cellfun(@(x) compareTimes(x), consSpks, 'UniformOutput', 0);
NSpkSpont = cellfun(@(x) sum(~x), evokSpksFlags);
trigs = trigs(trigs(:,2) <= round(consTime(2) * fs),:);
deltaTrigTime = sum([trigs(:,1); round(consTime(2) * fs)] - [0; trigs(:,2)])/fs;
firingRatePerCluster = NSpkSpont./deltaTrigTime;
end
