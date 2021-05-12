function [CV2, CVsqr, bIdx, Tot_CV2, Tot_CVsqr, Tot_bIdx, Tot_bTheta] =...
    compareFiringModes(condArray, spkSubs, fs, varargin)
%COMPAREFIRINGMODES computes the burst index per condition either in a
%block or trial by trial for the given spikes. The function requires the
%structure array conditions with 'name' and 'Triggers' (and 'Stimulus') in
%condArray, the spike times in indeces in spkSubs, and the sampling
%frequency in fs.
%   ? = compareFiringModes(condArray, spkSubs, fs)
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

getBurstsPerCluster = @(x, y) DiscreteWaveform.firstOfTrain(x{:}./fs, 10.^y);
[brsts, ~, ~, sps] = arrayfun(getBurstsPerCluster, spkSubs, gen_bTheta, fnOpts{:});
tnics = cellfun(@(x,y) xor(x,y), brsts, sps, fnOpts{:});

spkTms = cellfun(@(x) x./fs, spkSubs, fnOpts{:});
isiTms = cellfun(@(x) diff(x), spkTms, fnOpts{:});

if verbose
    % printing the threshold for non-bursty cells?
end

if blkFlag
    [CV2, CVsqr, bIdx] = compareConditionsInBlock();
else
    % not implemented yet
end

    function [CV2, CVsqr, bIdx] = compareConditionsInBlock()
        CV2 = zeros(Ncl, Ncond);
        CVsqr = CV2;
        bIdx = CV2;
        condSum = @(x,y) sum(x(y));
        for ccond = 1:Ncond
            % Getting the first trigger point and the last for the
            % considered condition in the array
            lims = [condArray(ccond).Triggers(1),...
                condArray(ccond).Triggers(end)];
            condBlockFlags = cellfun(@(x) x >= lims(1) & x <= lims(2),...
                spkSubs, fnOpts{:});
            brstsPerCluster = cellfun(condSum, brsts, condBlockFlags, fnOpts{:});
            tnicsPerCluster = cellfun(condSum, tnics, condBlockFlags, fnOpts{:});
            bIdx(:, ccond) = cellfun(@(x,y) x/(x+y),...
                brstsPerCluster, tnicsPerCluster);
            [CV2(:, ccond), CVsqr(:, ccond)] =...
                cellfun(@(x,y) getCVsfromISIs(x(y(1:numel(x)))), isiTms, condBlockFlags);
        end
    end

end