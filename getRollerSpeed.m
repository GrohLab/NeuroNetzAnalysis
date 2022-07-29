function [vf, rollTx] = getRollerSpeed(rollerPositions, fsRoll)
%GETROLLERSPEED computes the roller speed from the roller position table.
%   With rollerPositions and desired sampling frequency fsRoll, the funtion
%   interpolates asynchronous roller positions into a regularly sampled
%   signal. It differenciates and filters the results to get a smooth speed
%   signal.
%Emilio Isaias-Camacho@GrohLab 2021

us = 1e-6;
if isa(rollerPositions, "table")
    rt = [double(rollerPositions.RollerX), rollerPositions.RollerT];
    rt(circshift(diff(rt(:,2))==0, 1),:) = [];    
elseif isa(rollerPositions, "numeric")
    % Assuming the time resolution is micro seconds 
    rt = rollerPositions*diag([1,us]);
end
rt(isnan(rt(:,1)), :) = [];
rollTx = rt(1,2):1/fsRoll:rt(end,2); 
rxx = interp1(rt(:,2), rt(:,1), rollTx, "pchip");
% Filtering for 18 Hz
v = diff(rxx); [b, a] = butter(10, (2*18)/fsRoll, "low"); 
vf = filtfilt(b, a, v);
end