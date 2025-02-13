function [ rsq ] = goodnessFit2( pts, pts_hat, dim )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

SSEt = squeeze( sum( ( pts - pts_hat ).^2, dim ) );
SSTt = squeeze( sum( ( pts - mean( pts, dim ) ).^2, dim ) );
rsq = 1 - (SSEt./SSTt);

end

