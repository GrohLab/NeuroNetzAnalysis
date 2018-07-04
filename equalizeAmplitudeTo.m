function [outputSignal] = equalizeAmplitudeTo(referenceSignal,inputSignal)
%equalizeAmplitudeTo Summary of this function goes here
%   Detailed explanation goes here
if max(inputSignal) > 2*std(inputSignal)
    % If the mechanical TTL signal was active, normalize the amplitude
    % to the pressure signal.
    m = range(referenceSignal)/range(inputSignal);
    b = max(referenceSignal) - m*max(inputSignal);
    outputSignal = inputSignal*m + b;
else
    % If there was no mechanical stimulation, the signal would keep its
    % original amplitude but kept at the minimum values of the pressure
    % signal.
    outputSignal = (inputSignal - mean(inputSignal)) + min(referenceSignal);
end

