function [yhat, mdls] = fitSpline(x, y, n, winSz, ovrlap)
%FITSPLINE Fits an nth order piece-wise polynomial to the given signal y to
%the values inside the given window size winSz and overlapping the fittings
%ovrlap percentage (0-1)
%   [yhat, mdls] = fitSpline(x, y, n, winSz, ovrlap)
%       INPUTS:
%           x - vector containing the sample space of y (it frequently
%           represents time)
%           y - vector containing the values to be fitted
%           n - polynomial order
%           winSz - window size within the domain of x
%           ovrlap - percentage in which the windows would overlap to
%           average the fittings.
%       OUTPUTS:
%           yhat - spline estimation
%           mdls - n x R matrix containing the coefficients for each
%           fitting section, where n is the polynomial order and R is the
%           number of fittings done.
% Emilio Isa?as-Camacho @ GrohLab 2020
N = length(x);r = isrow(y); c = iscolumn(y);
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
dx = (winSz * (1 - ovrlap));
Nfits = range(x) / dx; dsubs = round(N/Nfits);
initSubs = ((0:ceil(Nfits)-1)*dsubs + 1)';
winSub = round((winSz * dsubs)/dx); ovrSub = winSub - dsubs;
mdls = zeros(n+1, ceil(Nfits)); yhat = nan(size(y));
dimy = [1,2]*[r;c];
ifVect = @(x) [x<N; x>=N];
%% Fitting a piecewise polynomial
for cft = 1:ceil(Nfits)
    subEnd = initSubs(cft)+winSub-1;
    csubs = initSubs(cft):[subEnd, N] * ifVect(subEnd);
    [mdl, yh] = fit_poly(x(csubs), y(csubs), n);
    
    yh = reshape(yh, 1*r + length(csubs)*c, 1*c + length(csubs)*r);
    mdls(:,cft) = mdl;
    ovrEnd = initSubs(cft)+ovrSub-1;
    coSubs = initSubs(cft):[ovrEnd, N] * ifVect(ovrEnd);
    ovrSub = [ovrSub, length(coSubs)] * ifVect(ovrEnd);
    avWeights = linspace(1, cft==1, ovrSub); 
    avWeights = reshape(avWeights, 1*r + ovrSub*c, 1*c + ovrSub*r);
    yh(1:ovrSub) =...
        sum(cat(dimy, yhat(coSubs).*avWeights,...
        yh(1:ovrSub).*flip(avWeights)), dimy, 'omitnan');
    yhat(csubs) = yh;
end
end
