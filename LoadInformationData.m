load('Data\120511c3.mat');
clearvars rep*
vt1 = V{1}(:,1)';
vt2 = V{2}(:,1)';
s = stim(:,1)';
sp1 = getSpikeTimes(vt1,sqrt(max(abs(vt1))),10);
sp2 = getSpikeTimes(vt2,sqrt(max(abs(vt2))),10);