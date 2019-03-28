%% PARAMETERS INITIALIZATION -- Modify this section for different experiments if necessary
% The timeLapse variable contains the time before the pulse in milliseconds
% (1 seconds = 1000 milliseconds & 1 millisecond = 0.001 second) and the
% time after the pulse in milliseconds. (30 ms before and 30 ms after)
timeLapse = [30,30];
% The laserPowers variable contains the series of laser intensities used to
% stimulate.
laserPowers = [0, 0.15, 2, 6, 10, 14, 20, 23];
% The laserTrials contains the number of pulses per laser intensity
% (Clustered into one figure)
laserTrials = 3;

%% Loading the recording
read_Intan_RHD2000_file
fprintf(1,'Thanks! Now filtering your data...\n')
fs = frequency_parameters.amplifier_sample_rate;
filteredSignals = iir50NotchFilter(amplifier_data',fs);
filteredSignals = iirSpikeFilter(filteredSignals,fs)';
[Nch, Ns] = size(filteredSignals);
fprintf(1,'Finished loading data\n')
%% Getting the condition trigger
dataCell = cell(32,1);
for ch = 1:32
    dataCell(ch) = {filteredSignals(ch,:)};
end
% These variables detect the rising and falling edges of the digital input
% pulse (laser).
stObj = StepWaveform(board_dig_in_data(4,:),fs,'On/Off','Laser');
% The lt variable contains a logical array size 2(on/off) x
% total_number_of_samples indicating the beginning and end of a pulse with
% a true value.
lt = stObj.Triggers;
% This contains the sample number at which the pulse rose. If you want to
% change for pulse fall, write as follows: ltOn = find(lt(:,2));
ltOn = find(lt(:,1));

[~,cSt] =...
    getStacks(false(1,length(amplifier_data(1,:))),...
    ltOn,'on',timeLapse * 1e-3,fs,fs,[],dataCell);
fprintf(1,'Got the condition triggers\n')

%% Plotting the continuous stack 
f = gobjects(8);
for laserPower = 1:numel(laserPowers)
    f(laserPower) =...
        figure('Name',sprintf('Power %.2f mW',laserPowers(laserPower)),...
        'Color',[1,1,1],'Renderer','painters','RendererMode','manual');
    for laserTrial = 1:laserTrials
        plotEEGchannels(...
            reshape(cSt(:,:,laserTrial + laserTrials*(laserPower - 1)),...
            Nch,sum(timeLapse*1e-3)*fs+1),1:Nch,sum(timeLapse*1e-3)+(1/fs),fs,0.5,f(laserPower));
        ax = f(laserPower).Children;
        ax.NextPlot = 'add';
    end
    [~,meanP] = plotEEGchannels(...
        reshape(mean(cSt(:,:,laserPower*3 - 2:laserPower*3),3),...
        Nch,sum(timeLapse*1e-3)*fs+1),1:32,sum(timeLapse*1e-3)+(1/fs),fs,0.5,f(laserPower));
    set(meanP,'Color',[0,0,0])
    savefig(f(laserPower),...
        fullfile(path,...
        sprintf('LaserStimulation_%.2fmW.fig',laserPowers(laserPower))))
end