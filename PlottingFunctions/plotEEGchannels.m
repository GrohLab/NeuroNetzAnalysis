function fhand = plotEEGchannels(EEG, labels, duration, fs, scale)
fhand = figure('Name','EEG');hold on
title('EEG channels')
[Nch, Ns] = size(EEG);
step = 5*std(EEG(:));
offset = step * Nch;
if fs > 1e9
    fprintf('Attempting a downsample of the signals only for displaying')
    fprintf(' purposes to 1 kHz\n')
    dFs = 1e3;
    DSfactor = round(fs/dFs);
    EEG = downsample(EEG',DSfactor)';
    Ns = size(EEG,2);
    dt = 1/dFs;
else
    dFs = fs;
    dt = 1/fs;
end
timeS = (0:Ns-1) * dt;
tick = zeros(1,Nch,'single');
tmSub = round(duration * dFs);
for c = 1:Nch
    tick(c) = offset;
    plot(timeS(1:tmSub),...
        scale*EEG(c,1:tmSub)+offset,'Color',repmat(0.8,1,3))
    offset = offset-step;
end
axis([0, timeS(tmSub),...
    tick(Nch)-2*step, tick(1)+2*step])
ylabel('Channels')
xlabel('Time [s]')
labels(strcmp(labels,'A1'))=[];
labels(strcmp(labels,'A2'))=[];
set(gca,'YTick',sort(tick),'YTickLabel',labels(Nch:-1:1))
end