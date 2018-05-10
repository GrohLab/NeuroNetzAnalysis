function [bIdx, tIdx, spIdx] = getInitialBurstSpike(spT,maxISI)
% GETINITIALBUSRTSPIKE gets the first spike time for each burst. The tonic
% spikes are unmodified according to the maximum inter-spiking-interval.
if ~isrow(spT)
    spT = spT';
end
delta_spT = cat(find(size(spT)~=1),inf,diff(spT));
% delta_spT = [inf,diff(spT)];
fsIdx = delta_spT >= maxISI;    % First spike index (burst-wise)
bsIdx = [~fsIdx(2:end) & fsIdx(1:end-1),fsIdx(end)];
tsIdx = ~bsIdx & fsIdx;
bIdx = spT(bsIdx);
tIdx = spT(tsIdx);
spIdx = spT(bsIdx | tsIdx);
end