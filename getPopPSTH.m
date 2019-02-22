function [psthStr] = getPopPSTH(dSt,condSt,configStruct)
%GETPOPPSTH returns the peri-stimulus triggered histogram for all the
%experiments collected in the population stack and conditions it according
%to both the configuration and conditioning structure.
%   Emilio Isaías-Camacho @GrohLab 2019

%% Initialization
% Gathering the necessary information to compute the PSTH combinations
% according to the user selection.
[Nv, Nts, Nap] = size(dSt.Stack);
fs = condSt.SamplingFrequencies(1); 
Ncv = size(condStr.CVDFlags,1); % Number of conditioning variables
finishedButton = false;
% To have a pure triggered PSTH (only taking random occurances of the
% ignored variables) we need to take those trials which were NOT excluded
% AND which had NO conditioning variables happening during the selected
% window of interest. De Morgan law: ~a & ~b = ~(a | b)
cleanIdx = ~(condSt.ExcludeFlags | any(condSt.CVDFlags,1));
cl_swaps = sum(cleanIdx);
fprintf('Percentage of ''clean'' trials: %.2f (%d/%d)\n',cl_swaps/Nap,...
    cl_swaps, Nap)
cleanPSTH = sum(dst.Stack(:,:,cleanIdx),3)./cl_swaps;

% If the user selected 1 conditioning variable, then there are no more
% possible combinations than the clean trigger and the trigger plus the
% conditioning variable.
if Ncv == 1
    uniquePSTH = sum(dst.Stack(:,:,~condSt.ExcludeFlags &...
        any(condSt.CVDFlags,1)),3);
else
    % Otherwise, the function will prompt the user to select a combination
    % of triggers in a list dialog for a later plotting.
    
    while ~finishedButton
        
    end
end

end

