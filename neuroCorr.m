function [corrStackCells] = neuroCorr(spks, tmReach, binSz, fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% getStacks(spT,alignP,ONOFF,timeSpan,fs,fsLFP,consEvents,...)
corrStackCells = cell(length(spks),1); Ncl = length(spks);
zeroCorr = ceil(tmReach * fs) + 1;
for ccl = 1:size(spks, 1)
    corrStackCells(ccl) = {squeeze(sum(getStacks(false, spks{ccl}, 'on',...
        [-tmReach, tmReach], fs, fs, spks(ccl+1:Ncl)),3))};
    if isempty(corrStackCells{ccl})
        Ne = size(ccl+1:Ncl,2); Nt = 2*ceil(tmReach*fs) + 1; 
        auxArray = zeros(Ne+2,Nt); m = memory; 
        Nspks = round((m.MaxPossibleArrayBytes * 0.9) / ((2+Ne) * Nt));
        Nspks_old = size(spks{ccl},1); spkProp = Nspks/Nspks_old;
        init = [(0:spkProp:1)';1]; edgeSubs = round(init * Nspks_old);
        fprintf(1, 'Separating the spikes in batches...\n')
        for csSet = 1:ceil(1/spkProp)
            fprintf(1,'Processing %d batch\n', csSet)
            ssSubs = (edgeSubs(csSet) + 1):edgeSubs(csSet + 1);
            auxArray = auxArray + squeeze(sum(getStacks(false,...
                spks{ccl}(ssSubs), 'on', [-tmReach, tmReach], fs, fs,...
                spks(ccl+1:Ncl)),3));
        end
        corrStackCells{ccl} = auxArray;
        fprintf(1, 'Done!\n')
    end
    corrStackCells{ccl}(2,:) = []; corrStackCells{ccl}(1,zeroCorr) = 0;
end
end
% corrStackCells = cellfun(@(x) squeeze(sum(getStacks(false, x, 'on',...
%     [-tmReach, tmReach], fs, fs, spks),3)), spks, 'UniformOutput', 0);