function [fanoPerCondition, tmAx] = computeFanoFactor(relSpkTms, binSize, varargin)
%COMPUTEFANOFACTOR estimates the intertrial variability for all given
%clusters.
%   Emilio Isa?as-Camacho @GrohLab

%% Parsing inputs
p = inputParser;

checkSpkTms = @(x) all([isstruct(x), isfield(x, {'name','SpikeTimes'})]);
checkBinSize = @(x) all([isnumeric(x), x > 0]);

defTmWin = [-0.1, 0.1];
checkTmWin = @(x) all([numel(x) == 2, x(1) < x(2)]);

defKernelLength = 10e-3;
checkKernLen = @(x) all([x > binSize, isnumeric(x), x > 0, numel(x) == 1]);

defKernel = @rectwin;
checkKern = @(x) all([isa(x, 'function_handle'),...
    strfind(func2str(x),'win') | strcmp(func2str(x),'triang') |... 
    strfind(func2str(x),'blackman') | strcmp(func2str(x),'hamming') |...
    strcmp(func2str(x),'hann') | strcmp(func2str(x),'kaiser') | ...
    strcmp(func2str(x),'bartlett')]);

defKernSpecial = {};

defAdjustTmWinFlag = false;
checkAdjustFlag = @(x) all([isnumeric(x) | islogical(x), numel(x) == 1]);

defVerbose = false;
checkVerbose = @(x) all([isnumeric(x) | islogical(x), numel(x) == 1]);

p.addRequired('relSpkTms', checkSpkTms);
p.addRequired('binSize', checkBinSize);

p.addParameter('timeWindow', defTmWin, checkTmWin);
p.addParameter('kernelLength', defKernelLength, checkKernLen);
p.addParameter('kernel', defKernel, checkKern);
p.addParameter('kernelParameters', defKernSpecial, @iscell);
p.addParameter('adjustTime', defAdjustTmWinFlag, checkAdjustFlag);

p.addOptional('verbose', defVerbose, checkVerbose);

p.KeepUnmatched = true;

p.parse(relSpkTms, binSize, varargin{:})

relSpkTms = p.Results.relSpkTms;
binSize = p.Results.binSize;
tmWin = p.Results.timeWindow;
kernLen = p.Results.kernelLength;
kernfun = p.Results.kernel;
kernParams = p.Results.kernelParameters;
adjustFlag = p.Results.adjustTime;
verb = p.Results.verbose;

fnOpts = {'UniformOutput', false};

%% Computation of the fano factor
% FF = var(#spikes)/mean(#spikes)

Ncond = size(relSpkTms,2);

%     function [dataTmWin] = getTimeLimitsFromData(relSpkTms)
%         auxIdx = 1;
%         dataTmWin = zeros(2,1);
%         for func = {@min, @max}
%             auxTm = arrayfun(@(x) cellfun(func{:},...
%                 relSpkTms(x).SpikeTimes, fnOpts{:}), (1:Ncond)', fnOpts{:});
%             while iscell(auxTm)
%                 auxTm = cat(2, auxTm{:});
%             end
%             dataTmWin(auxIdx) = func{:}(auxTm); auxIdx = auxIdx + 1;
%         end
%     end

% dataTmWin = getTimeLimitsFromData(relSpkTms);

histOpts = {'BinWidth', binSize, 'BinLimits', tmWin};
getFanoFactor = @(x) var(x, 'omitnan')./mean(x, 'omitnan');
Nclust = size(relSpkTms(1).SpikeTimes, 1);
fanoPerCondition = cell(Ncond,1);
Nbins = range(tmWin)/binSize;
winLen10ms = 2^nextpow2(round(kernLen/binSize));
winLen = winLen10ms;
if winLen10ms > Nbins
    winLen = Nbins/5;
end
if winLen10ms < 4
    winLen = 4;
end
try
    kern = kernfun(winLen, kernParams{:});
catch
    if verb
        fprintf(1, '%s doesn''t accept more parameters!\n',...
            func2str(kernfun))
    end
    kern = kernfun(winLen);
end

for ccond = 1:Ncond
    conditionName = relSpkTms(ccond).name;
    if verb
        fprintf(1, 'Computing fano factor for %s (%d/%d)\n', conditionName,...
            ccond, Ncond)
    end
    spkTmsCell = relSpkTms(ccond).SpikeTimes;
    [trialPsth, tmAx] = arrayfun(@(x) cellfun(@(y) histcounts(y, histOpts{:}),...
        spkTmsCell(x,:), fnOpts{:}), (1:Nclust)', fnOpts{:});
    trialPsth = cellfun(@(x) cat(1, x{:}), trialPsth, fnOpts{:});
    trialPsth = cellfun(@(x) conv2(x', kern, 'same')', trialPsth, fnOpts{:});
    trialPsth = cellfun(getFanoFactor, trialPsth, fnOpts{:});
    fanoPerCondition{ccond} = cat(1, trialPsth{:});
end
while iscell(tmAx)
    tmAx = tmAx{1};
end
tmAx = arrayfun(@(x) mean(tmAx(x:x+1)), (1:numel(tmAx)-1)');
end

