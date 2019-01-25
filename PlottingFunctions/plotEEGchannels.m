function fhand = plotEEGchannels(EEG, labels, duration, fs, scale)
fhand = figure('Name','EEG');hold on
title('EEG channels')
ch=size(EEG,1);
step = 5*std(EEG(:));
offset = step * ch;
N = size(EEG,2);
timeS = (0:N-1)*(1/fs);
tick = zeros(1,ch,'single');
for c = 1:ch
    tick(c) = offset;
    plot(timeS(1:1+duration*fs),...
        scale*EEG(c,1:1+duration*fs)+offset,'b')
    offset = offset-step;
end
axis([0, timeS(1+duration*fs),...
    tick(ch)-2*step, tick(1)+2*step])
ylabel('Channels')
xlabel('Time [s]')
labels(strcmp(labels,'A1'))=[];
labels(strcmp(labels,'A2'))=[];
set(gca,'YTick',sort(tick),'YTickLabel',labels(ch:-1:1))
end