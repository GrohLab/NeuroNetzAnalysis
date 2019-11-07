function [H, P] = statTests(obs1, obs2, varargin)
%STATTESTS performs the requested test between the two given observations
%sets. The possible observations are the following: exact binomial,
%McNemar, Chi^2, kstest.
%
%   [H, P] = statTests(obs1, obs2, testType)
%
%       INPUT
%           - obs1 -- vector or matrix containing the first set of
%           observations
%           - obs2 -- vector or matrix containing the second set of
%           observations
%           - testType -- character vector describing the test to be
%           performed ('binomial', 'mcnemar', 'chi2', 'kstest')
%       OUTPUT
%           - H -- logical vector or matrix indicating whether the compared
%           observations reject the null hypothesis; if the P value is
%           under 0.05
%           - P -- vector or matrix containing the significance values.



end

