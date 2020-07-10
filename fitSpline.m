function [yhat, mdls] = fitSpline(x, y, n, winSz, ovrlap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = length(x);
if winSz >= max(x)
    fprintf(1, 'Window is bigger than the signal!\n')
    fprintf(1, 'Setting the window size to the complete signal length (%d)',...
        max(x))
    winSz = max(x); ovrlap = 0;
end
if n >= 1 
    fprintf(1, 'Fitting a polynomial of %d degree...\n', round(n))
else
    warning('''n'' should be a positive integer')
    fprintf(1,'Setting ''n'' to 1...\n')
    n = 1;
end
if ovrlap >= 1 || ovrlap < 0
    warning('The overlap should be a number greater or equal to zero and smaller than one!')
    fprintf(1,'Setting the overlap to 35%%\n')
    ovrlap = 0.35;
end
Nfits = ceil(max(x) / (winSz * (1 - ovrlap)));
initSubs = (0:Nfits-1)*(winSz * (1 - ovrlap));
%% Fitting a piecewise polynomial
for cft = 1:Nfits
    
end
end
