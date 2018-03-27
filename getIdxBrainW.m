%function [ alpha_idxs, beta_idxs, delta_idxs, theta_idxs, total_idxs] =...
function [ idxsM ] = getIdxBrainW(freq_axis,brain_bands)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if isempty(brain_bands)
    brain_bands = [7,13;13.1,30.0;0.5,4;4.1,6.9;0.5,35.0];
end
% alpha_idxs = freq_axis >= brain_bands(1,1) & freq_axis <= brain_bands(1,2);
% beta_idxs  = freq_axis >= brain_bands(2,1) & freq_axis <= brain_bands(2,2);
% delta_idxs = freq_axis >= brain_bands(3,1) & freq_axis <= brain_bands(3,2);
% theta_idxs = freq_axis >= brain_bands(4,1) & freq_axis <= brain_bands(4,2);
% total_idxs = freq_axis >= brain_bands(5,1) & freq_axis <= brain_bands(5,2);
[Nb, ~] = size(brain_bands);
for bw = 1:Nb
    idxsM(bw,1) = find(freq_axis >= brain_bands(bw,1),1);
    idxsM(bw,2) = find(freq_axis <= brain_bands(bw,2),1,'last');
end
end