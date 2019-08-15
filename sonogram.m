function [sonoStruct] =...
    sonogram(signal,sampling_frequency, overlap, window_duration)
%SONOGRAM Returns the time-frequency components for the given signal.
%   This function accepts 4 arguments: signal, which is the time series to
%   be transformed; sampling_frequency of the given signal; window_size,
%   which is the period in seconds for the Fourier transform to be applied;
%   and overlap, which is the percentage for the moving windows to retake
%   into the FFT from the previous one.
%   Emilio Isaias-Camacho @ GrohLab 2018
N = length(signal);
if nargin == 3
    Nw = 2^log2(sampling_frequency);
elseif nargin == 4
    Nw = round(window_duration * sampling_frequency);
end
fobj = FourierSpectrum(zeros(1,Nw,'single'),sampling_frequency);
Nzp = fobj.N;
NonOverWin = round(Nw * (1-overlap));
lenImg = round(N/NonOverWin);
SpectrumImage = zeros((Nzp/2),lenImg,'single');
% MagnitudeImage = zeros((Nzp/2),lenImg);
% AngleImage = MagnitudeImage;
currentWindow = 1;
cw = 1;
h = waitbar(0,'Computing the Sonogram...');
while currentWindow + Nw <= N
    sigSeg = signal(currentWindow:currentWindow+Nw-1);
    fobj = FourierSpectrum(sigSeg,sampling_frequency);
    SpectrumImage(:,cw) = fobj.getHalFourier;
    currentWindow = currentWindow + NonOverWin;
    cw = cw + 1;
    waitbar(currentWindow/N)
end
close(h)
tx = 0:NonOverWin/sampling_frequency:(NonOverWin*(lenImg-1))/sampling_frequency;
fx = 0:sampling_frequency/(Nzp):(sampling_frequency*(Nzp-2))/(2*Nzp);
sonoStruct = struct('SpectrumImage',SpectrumImage,...
    'TimeAxis',tx,'FrequencyAxis',fx);
end

