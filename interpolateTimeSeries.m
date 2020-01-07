function [interpVals] = interpolateTimeSeries(timeVector, yValues, desTime)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Ny = numel(desTime);
interpVals = zeros(Ny,1);
for cdt = 1:numel(desTime)
    if desTime(cdt) > max(timeVector) || desTime(cdt) < min(timeVector)
        interpVals(cdt) = NaN;
        continue
    end
    [~, clsts] = min(abs(timeVector - desTime(cdt)));
    lnTime = timeVector(clsts:clsts+1);
    lnYval = yValues(clsts:clsts+1);
    if timeVector(clsts) > desTime(cdt)
        lnTime = timeVector(clsts-1:clsts);
        lnYval = yValues(clsts-1:clsts);
    end
    mdl = fit_poly(lnTime, lnYval, 1);
    interpVals(cdt) = mdl(1)*desTime(cdt) + mdl(2);
end
end

