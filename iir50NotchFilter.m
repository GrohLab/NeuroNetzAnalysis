function [csig] = iir50NotchFilter(signal,sampling_frequency)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% w0 = cut_frequency / (sampling_frequency/2)

w0 = 100/sampling_frequency;
% bw is the bandwidth of the notch using a 'Q' factor.
bw = w0/35;
[b50, a50] = iirnotch(w0,bw);
fprintf('Filtering the signal for 50 Hz...')
csig = filtfilt(b50,a50,signal);
fprintf(' done!\n');
end