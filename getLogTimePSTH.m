function [PSTHstruct] = getLogTimePSTH(relSpkTms, responseFlags, varargin)
%LOGTIMEPSTH creates the PSTH in log time in between a given time window
%using the relative spike times structure and boolean flags to identify
%responsive clusters.
%   [outputArg1] = logTimePSTH(relSpkTms, tmWin, Name, Value)
%   INPUTS:
%       relSpkTms - structure array with 'name' and 'SpikeTimes' as fields.
%                   Produced in DE_Jittering and poolExperiments scripts.
%       responseFlags - NTclx1 logical vector identifying which clusters
%                       are responding to the stimulus.
%       NAME  -- VALUE
%       tmWin - 2 element vector containing time window edges for
%               considering the log-PSTH. Warning! The minimal value cannot
%               be lower or equal than zero. Log(0) = -inf. Default 1 to 50
%               ms.
%       fs - sampling frequency. Default 1/33.3e-6 s ~30030.03 Hz
%       Nbin - number of logarithmic bins to divide the time window.
%              Default 64.
%       *verbose - [optional] logical flag to activate (1) or deactive (0 -
%                  default)
%       *Offset - [optional] shift in time in seconds for the spikes.
%   OUTPUTS:
%       PSTHstruct - structure with 'LogPSTH' which is a  Cl x Nbin x Cd
%       matrix containing the activity for Cl clusters, Nbin bins, and Cd
%       conditions; 'Log10TimeAxis' which is a 1 x Nbin vector containing
%       the exponent base 10; 'ConditionNames' requires no explanation;
%       'DeltaLogStep' which is the step in log domain; and 'TimeAxis'
%       which is a 1 x Nbin vector containing the bin center in seconds.
%
%Emilio Isaias-Camacho @ GrohLab 2021

%% Checking the input arguments
p = inputParser;
% Required arguments
checkSpkStruct = @(x) any([isstruct(x), isfield(x,{'name','SpikeTimes'})]);
checkResponseFlag = @(x) any([size(x,1) == size(relSpkTms(1).SpikeTimes,1),...
    isnumeric(x) | islogical(x)]);
% Required parameters
defTmWin = [0.001, 0.05];
checkTmWin = @(x) all([x > 0, diff(x) > 0, numel(x) == 2]);

defFs = 1/33.3e-6;
checkFs = @(x) all([isnumeric(x), x>0]);

defNbin = 64;
checkNbin = @(x) all([isnumeric(x), x>1]);

% Optional parameter
defVerbose = false;
checkVerbose = @(x) all([numel(x) == 1, isnumeric(x) | islogical(x)]);

defOffset = 0;
checkOffset = @(x) all([isnumeric(x), ~isnan(x), ~isinf(x)]);

p.addRequired('relSpkTms', checkSpkStruct);
p.addRequired('responseFlags', checkResponseFlag);
p.addParameter('tmWin', defTmWin, checkTmWin);
p.addParameter('fs', defFs, checkFs);
p.addParameter('Nbin', defNbin, checkNbin);
p.addOptional('verbose', defVerbose, checkVerbose);
p.addOptional('Offset', defOffset, checkOffset);

p.KeepUnmatched = true;

p.parse(relSpkTms, responseFlags, varargin{:});

relSpkTms = p.Results.relSpkTms;
responseFlags = p.Results.responseFlags;
tmWin = p.Results.tmWin;
fs = p.Results.fs;
Nbin = p.Results.Nbin;
verbose = p.Results.verbose;
ofst = p.Results.Offset;

if verbose
    fprintf(1, 'logPSTH from %.3f to %.3f ms with %d bins\n',...
        tmWin*1e3, Nbin)
end
%% Auxiliary variables
fnOpts = {'UniformOutput', false};
hstOpts = {'Normalization','probability'};
Ncond = numel(relSpkTms);
conditionNames = arrayfun(@(x) x.name, relSpkTms, fnOpts{:});
[binCenters, binEdges, ~, deltaLogT] = prepareLogBinEdges(tmWin, Nbin);
% Matrix holding the logbin counts.
logPSTH = zeros(sum(responseFlags), Nbin, Ncond);
if verbose
    fprintf(1, '%d responsive clusters out of %d in structure\n',...
        sum(responseFlags), size(relSpkTms(1).SpikeTimes,1))
end

%% Main loop
if verbose
    fprintf(1, 'Starting main loop... ')
end
for ccond = 1:Ncond
    spkTms = arrayfun(@(x) cat(2, relSpkTms(ccond).SpikeTimes{x,:})+ofst,...
        (1:sum(responseFlags))', fnOpts{:});
    spkTms = cellfun(@(x) x(x > tmWin(1) & x <= tmWin(2)), spkTms,...
        fnOpts{:});
    binCount = cellfun(@(x) histcounts(log10(x), binEdges, hstOpts{:}),...
        spkTms, fnOpts{:}); binCount = cat(1, binCount{:});
    logPSTH(:,:,ccond) = binCount;
    if verbose
        fprintf(1, '%s processed, ', conditionNames{ccond})
    end
end
PSTHstruct = struct('LogPSTH',logPSTH, 'Log10TimeAxis', binCenters,...
    'TimeAxis', 10.^binCenters, 'ConditionNames', string(conditionNames),...
    'DeltaLogStep', deltaLogT);
if verbose
    fprintf(1, ' done!\n')
end
end