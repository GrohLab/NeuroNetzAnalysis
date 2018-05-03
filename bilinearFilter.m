function [f_new] = bilinearFilter(wind_size, signal)
% BILINEARFILTER performs a filtering operation which preserves the edges
% information and reduces the noise in flat areas.
%   [smth, flt] = bilinearFilter(window_size, signal)

cw = floor(wind_size/2) +1;
f_new = zeros(size(signal));
s = @(f1,f2) (exp((-(abs(repmat(f1,size(f2))-f2)/std(f2))/2).^2));
N = length(signal);
winHalf = floor(wind_size/2);
h = waitbar(0,'Bilinear filter...');
while cw+3 + winHalf <= N
    similarity1 = s(signal(cw),signal(cw-winHalf:cw+winHalf));
    similarity2 = s(signal(cw+1),signal(cw+1-winHalf:cw+1+winHalf));
    similarity3 = s(signal(cw+2),signal(cw+2-winHalf:cw+2+winHalf));
    similarity4 = s(signal(cw+3),signal(cw+3-winHalf:cw+3+winHalf));
    closeness1 = s(cw,cw-winHalf:cw+winHalf);
    closeness2 = s(cw+1,cw+1-winHalf:cw+1+winHalf);
    closeness3 = s(cw+2,cw+2-winHalf:cw+2+winHalf);
    closeness4 = s(cw+3,cw+3-winHalf:cw+3+winHalf);
    f_new(cw) = sum(similarity1.*signal(cw-winHalf:cw+winHalf).*...
        closeness1)/(similarity1*closeness1');
    f_new(cw+1) = sum(similarity2.*signal(cw+1-winHalf:cw+1+winHalf).*...
        closeness2)/(similarity2*closeness2');
    f_new(cw+2) = sum(similarity3.*signal(cw+2-winHalf:cw+2+winHalf).*...
        closeness3)/(similarity3*closeness3');
    f_new(cw+3) = sum(similarity4.*signal(cw+3-winHalf:cw+3+winHalf).*...
        closeness4)/(similarity4*closeness4');
    cw = cw + 4;
    waitbar(cw/(N-winHalf))
end
close(h)
% s = zeros(1,wind_size);
% c = s;
% flt = zeros(size(signal));
% smth = flt;
% sig = std(signal);
% for ii = 1:length(signal)-wind_size+1
%     auxseg = signal(ii:ii+wind_size-1);
%     for jj = 1:length(auxseg)
%         s(jj) = exp(-0.5* (...
%             norm(auxseg(jj) - auxseg(round(wind_size/2)))/sig)^2);
%         c(jj) = exp(-0.5* (...
%             norm(jj - round(wind_size/2))/round(wind_size/2))^2);
%     end
%     s = s/norm(s);
%     c = c/norm(c);
%     smth(ii+round(wind_size/2)-1) = (s.*c)*auxseg;
%     flt(ii+round(wind_size/2)-1) = s*c';
% end
% 
% end