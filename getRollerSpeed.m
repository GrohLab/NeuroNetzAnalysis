function [vf, rollTx] = getRollerSpeed(rollerPositions, fsRoll)
%GETROLLERSPEED computes the roller speed from the roller position table.
%   With rollerPositions and desired sampling frequency fsRoll, the funtion
%   interpolates asynchronous roller positions into a regularly sampled
%   signal. It differenciates and filters the results to get a smooth speed
%   signal.
%Emilio Isaias-Camacho@GrohLab 2021

% Assuming the time resolution is micro seconds 
rollTx = rollerPositions(1,2)/1e6:1/fsRoll:rollerPositions(end,2)/1e6;
rxx = interp1(rollerPositions(:,2)/1e6, rollerPositions(:,1), rollTx, "pchip");
% Filtering for 18 Hz
v = diff(rxx); [b, a] = butter(10, (2*18)/fsRoll, "low"); 
vf = filtfilt(b, a, v);
end