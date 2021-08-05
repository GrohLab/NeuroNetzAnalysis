function [CV2, CVsqr, bIdx, gen_bTheta, brstPack,...
    Tot_CV2, Tot_CVsqr, Tot_bIdx, Tot_bTheta, bCounts] =...
    compareFiringModes(condArray, spkSubs, fs, varargin)
%COMPAREFIRINGMODES computes the burst index per condition either in a
%block or trial by trial for the given spikes. The function requires the
%structure array conditions with 'name' and 'Triggers' (and 'Stimulus') in
%condArray, the spike times in indeces in spkSubs, and the sampling
%frequency in fs. It returns the coefficient of variation, CV^2, burst
%index, and the used threshold for each cluster in each condition (per
%block or trial) and regardless of the conditions.
%   [CV2, CVsqr, bIdx, bTh] = compareFiringModes(condArray, spkSubs, fs)
%   compareFiringModes(condArray, spkSubs, fs, Name-Value)
%
%   INPUTS:
%       condArray - vector array with C elements for the considered
%                   conditions to look at.
%       spkSubs - vector cell array with N elements each containing the
%                 spike subscripts representing the time of occurrance.
%       fs - scalar indicating the sampling frequency for the given spikes
%   NAME-VALUE:
%       *inBlock - [optional] boolean flag indicating if the conditions
%                   were in a block or in as independent and interleaved
%                   trials. (default: block).
%       *verbose - [optional] boolean flag indicating if the function
%                   should display messages to the user
%   OUTPUTS:
%       CV2, CVsqr, bIdx, bTh - NxC matrices containing 3 measurements for
%                               N neurons and C conditions, and the
%                               considered threshold: coefficient of
%                               variation, CV^2, burst index (burst
%                               proportion), and maximum inter-spike
%                               interval to consider a putative burst.
%       totCV2, totCVsqr, totbIdx - Nx1 matrices containing  measurements
%                                   for N neurons in general.

%Emilio Isaias-Camacho @GrohLab 2021
%% Parsing inupts
p = inputParser;

checkCondArray = @(x) all([isstruct(x), isfield(x,{'name','Triggers'})]);

checkFs = @(x) all([isnumeric(x), numel(x) == 1, x > 0]);

defBlkFlg = true;
checkInBlock = @(x) all([isnumeric(x) | islogical(x), numel(x) == 1]);

defVerbose = false;
checkVerbose = checkInBlock;

p.addRequired('condArray', checkCondArray);
p.addRequired('spkSubs', @iscell);
p.addRequired('fs', checkFs);

p.addParameter('inBlock', defBlkFlg, checkInBlock);
p.addOptional('verbose', defVerbose, checkVerbose);

p.parse(condArray, spkSubs, fs, varargin{:});

condArray = p.Results.condArray;
spkSubs = p.Results.spkSubs;
fs = p.Results.fs;
blkFlag = p.Results.inBlock;
verbose = p.Results.verbose;

%% Compute the total mode for all spikes
fnOpts = {'UniformOutput', false};

[Tot_CV2, Tot_CVsqr, Tot_bIdx, Tot_bTheta] = getBurstInfo(spkSubs, fs);
genericTheta = nanmean(Tot_bTheta(Tot_bTheta < -2)); gen_bTheta = Tot_bTheta;
gen_bTheta(Tot_bTheta > log10(2)-2) = genericTheta;

Ncond = size(condArray(:),1);
Ncl = size(spkSubs(:),1);
% First spike in a burst
getBurstsPerCluster = @(x, y) DiscreteWaveform.firstOfTrain(x{:}./fs, 10.^y);
[brsts, ~, ~, sps] = arrayfun(getBurstsPerCluster, spkSubs, gen_bTheta, fnOpts{:});
% Last spike in a burst
getLstBstPerCluster = @(x, y) DiscreteWaveform.lastOfTrain(x{:}./fs, 10.^y);
lstbst = arrayfun(getLstBstPerCluster, spkSubs, gen_bTheta, fnOpts{:});

% Counting spikes per burst
fstSpk = cellfun(@find, brsts, fnOpts{:});
lstSpk = cellfun(@find, lstbst, fnOpts{:});
bCounts = brsts;
for ccl = 1:Ncl
    bCounts{ccl}(brsts{ccl}) = (lstSpk{ccl} - fstSpk{ccl}) + 1;
end
% Tonic spikes (or with greater ISI than the threshold)
tnics = cellfun(@(x,y) xor(x,y), brsts, sps, fnOpts{:});

spkTms = cellfun(@(x) x./fs, spkSubs, fnOpts{:});
isiTms = cellfun(@(x) diff(x), spkTms, fnOpts{:});

if verbose
    % printing the threshold for non-bursty cells?
end

if blkFlag
    [CV2, CVsqr, bIdx, brstPack] = compareConditionsInBlock();
else
    % not implemented yet
end

    function [CV2, CVsqr, bIdx, brstPack] = compareConditionsInBlock()
        CV2 = zeros(Ncl, Ncond);
        CVsqr = CV2;
        bIdx = CV2;
        condSum = @(x,y) sum(x(y));
        brstPack = cell(Ncl, Ncond);
        sliceLogicalTrace = @(x,y) x(y);
        for ccond = 1:Ncond
            % Getting the first trigger point and the last for the
            % considered condition in the array
            lims = [condArray(ccond).Triggers(1),...
                condArray(ccond).Triggers(end)];
            condBlockFlags = cellfun(@(x) x >= lims(1) & x <= lims(2),...
                spkSubs, fnOpts{:});
            brstsPerCluster = cellfun(condSum, brsts, condBlockFlags, fnOpts{:});
            tnicsPerCluster = cellfun(condSum, tnics, condBlockFlags, fnOpts{:});
            brstPack(:,ccond) = cellfun(sliceLogicalTrace, bCounts,...
                condBlockFlags, fnOpts{:});
            bIdx(:, ccond) = cellfun(@(x,y) x/(x+y),...
                brstsPerCluster, tnicsPerCluster);
            [CV2(:, ccond), CVsqr(:, ccond)] =...
                cellfun(@(x,y) getCVsfromISIs(x(y(1:numel(x)))), isiTms, condBlockFlags);
        end
        brstPack = cellfun(@(x) x(x~=0), brstPack, fnOpts{:});
    end

end