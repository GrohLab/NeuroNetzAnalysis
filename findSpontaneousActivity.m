function [spontaneouStack] =...
    findSpontaneousActivity(spT, discreteTraces, timeLapse, fs, per100)
%FINDSPONTANEOUSACTIVITY returns the activity in which there is no active
%conditioning variables in the experiment. It is structured in a stacked
%form and a PSTH can be easily computed.
%   More information will be given in the future
[Ne, Ns] = size(discreteTraces);
spontaneousTrace = sum(discreteTraces,1) < 1;
sponSteps = StepWaveform(spontaneousTrace, fs);
sponBouts = sponSteps.Triggers;
sponIdx = diff(sponBouts,1,2)/fs * (1-per100) >= sum(timeLapse);
end

