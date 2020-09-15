function [ptsW, Wd] = whitenPoints(pts)
%WHITENPOINTS computes the whitening matrix of and transforms the given
%points without their mean.

% Remove the mean
pts = pts - mean(pts);
% Compute the covariance matrix
sig = cov(pts);
% Compute the whitening matrix
[E, D] = eig(sig);
D = pinv(real(sqrt(D))); % 1/sqrt(D)
Wd = E*D;
% Apply the transformation
ptsW = pts*Wd;
end