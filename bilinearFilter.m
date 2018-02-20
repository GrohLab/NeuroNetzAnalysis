function [smth, flt] = bilinearFilter(wind_size, signal)
% BILINEARFILTER performs a filtering operation which preserves the edges
% information and reduces the noise in flat areas.
%   [smth, flt] = bilinearFilter(window_size, signal)
s = zeros(1,wind_size);
c = s;
flt = zeros(size(signal));
smth = flt;
sig = std(signal);
for ii = 1:length(signal)-wind_size+1
    auxseg = signal(ii:ii+wind_size-1);
    for jj = 1:length(auxseg)
        s(jj) = exp(-0.5* (...
            norm(auxseg(jj) - auxseg(round(wind_size/2)))/sig)^2);
        c(jj) = exp(-0.5* (...
            norm(jj - round(wind_size/2))/round(wind_size/2))^2);
    end
    s = s/norm(s);
    c = c/norm(c);
    smth(ii+round(wind_size/2)-1) = (s.*c)*auxseg;
    flt(ii+round(wind_size/2)-1) = s*c';
end

end