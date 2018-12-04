function [outputSpectrum,desFx] = interpolateSpectrum(inputSpectrum,halfFlag,inputFs,desiredFs,desN)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N = length(inputSpectrum);
if halfFlag
    deltaOmega = inputFs/(2*length(inputSpectrum));
    desDeltaOmega = desiredFs/(2*desN);
    fx = (0:N-1)*deltaOmega;
    desFx = (0:N-1) * desDeltaOmega;
else
    deltaOmega = inputFs/length(inputSpectrum);
    desDeltaOmega = desiredFs/desN;
    fx = (0:N-1)*deltaOmega - inputFs/2;
    desFx = (0:N-1) * desDeltaOmega - desiredFs/2;
end

fprintf('The resolution of the interpolated spectrum is %.5f Hz/div.\n',...
    desDeltaOmega)
outputSpectrum = zeros(1,desN);
h = waitbar(1,'Interpolating spectrum');
for cs = 2:desN
    x = desFx(cs);
    [~,x1sub]=min(distmatrix(fx',x));
    if fx(x1sub) > x
        x2sub = x1sub;
        x1sub = x1sub - 1;
    elseif x1sub < N
        x2sub = x1sub + 1;
    else
        x2sub = x1sub;
        deltaOmega = 1;
    end
    x1 = fx(x1sub);x2 = fx(x2sub);
    y1 = inputSpectrum(x1sub);y2 = inputSpectrum(x2sub);
    outputSpectrum(cs) = (y1*distmatrix(x,x2) + y2*distmatrix(x,x1))/...
        deltaOmega;
    waitbar(cs/desN)
end
outputSpectrum(1) = inputSpectrum(1);
close(h)
end

