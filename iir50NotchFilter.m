function [csig] = iir50NotchFilter(signal,sampling_frequency)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
w0 = 100/sampling_frequency;
bw = w0/35;
[b50, a50] = iirnotch(w0,bw);
fprintf('Filtering the signal for 50 Hz...')
csig = filtfilt(b50,a50,signal);
fprintf(' done!\n');
end

