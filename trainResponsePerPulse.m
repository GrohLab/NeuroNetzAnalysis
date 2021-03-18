function [ptsOn,ptsOf,sOn,sOf] =...
    trainResponsePerPulse(PSTH, txpsth, freq, dur, Npulse, interstngWin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dt = 1/freq;
onst = (0:Npulse-1)'*dt;
ofst = (0:Npulse-1)'*dt + dur;

onrpWins = [onst, onst] + interstngWin;
ofrpWins = [ofst, ofst] + interstngWin;
onrpIdx = txpsth >= onrpWins(:,1) & txpsth <= onrpWins(:,2);
ofrpIdx = txpsth >= ofrpWins(:,1) & txpsth <= ofrpWins(:,2);
ptsOn = zeros(size(onrpIdx,1),size(PSTH,2),2); % time, magnitude
ptsOf = ptsOn;
sOn = size(Npulse, size(PSTH,2)); sOf = sOn;
for ccond = 1:size(PSTH,2)
    for crw = 1:size(onst,1)
        [mg, tmSub] = max(PSTH(onrpIdx(crw,:),ccond));
        tmWinSub = find(onrpIdx(crw,:));
        ptsOn(crw, ccond, 1) = txpsth(tmWinSub(tmSub));
        ptsOn(crw, ccond, 2) = mg;
        [mg, tmSub] = max(PSTH(ofrpIdx(crw,:),ccond));
        tmWinSub = find(ofrpIdx(crw,:));
        ptsOf(crw, ccond, 1) = txpsth(tmWinSub(tmSub));
        ptsOf(crw, ccond, 2) = mg;
        sOn(crw) = sum(PSTH(onrpIdx(crw,:),ccond),1);
        sOf(crw) = sum(PSTH(ofrpIdx(crw,:),ccond),1);
    end
end
end

