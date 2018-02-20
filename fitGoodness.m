function [ r2 ] = fitGoodness( poly_fit,y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
SS=sumsqr(y-mean(y));
r2 = 1 - (sumsqr(y - poly_fit)/SS);
end

