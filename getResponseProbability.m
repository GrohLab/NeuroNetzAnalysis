function [p_f, pEvok, pSpon, fishStats] = ...
    getResponseProbability(relSpkTms,cStruct,varargin)
%GETRESPONSEPROBABILITY The trial is classified as responsive (true) when
%the spiking activity within the response window exceeds the spiking
%activity in the spontaneous window
%   [p_f] = getResponseProbability(relSpkTms,cStruct)
%   [p_f] = getResponseProbability(relSpkTms,cStruct,varargin)
%   INPUTS:
%       relSpkTms - 1xC structure array with 'name' and 'SpikeTimes' as
%                   fields. Produced in DE_Jittering and poolExperiments
%                   scripts.
%       cStruct - structure containing the information of the 
%   NAME-VALUE:
%       'sponWin' - 1x2 array representing the time when the neuron is
%                   expected to respond. Default is -flip(respWin)
%       'timeLimits' - 1x2 array representing the time for limiting the
%                      time histograms. Default is -0.1 to 0.1 seconds.
%       'binWidth' - scalar value specifying the bin size for the time
%                    histograms in seconds. Default is 1 ms (1e-3 s).
%       'verbose' - logical value for printing use feedback. Default is
%                   false.
%   OUTPUT:
%       p_f - NxC array with the probability of response for each neuron N
%             and condition C. 
%       pSpon - NxC array with the probability of spiking occurrances
%               within the spontaneous window
%       pEvok - NxC array with the probability of spiking occurrances
%               within the response window
%Emilio Isaias-Camacho @ GrohLab 2021

%% Parse inputs

pfs = inputParser;
% Required arguments
checkSpkStruct = @(x) all([isstruct(x), isfield(x,{'name','SpikeTimes'})]);

% Optional parameter
%   Verbose
defVerbose = false;
checkVerbose = @(x) all([numel(x) == 1, isnumeric(x) | islogical(x)]);

% Parsing
pfs.addRequired('relSpkTms', checkSpkStruct);
pfs.addRequired('cStruct', @checkConfigStruct);
pfs.addOptional('verbose', defVerbose, checkVerbose); 
pfs.KeepUnmatched = true;
pfs.parse(relSpkTms, cStruct, varargin{:});

% Results
relSpkTms = pfs.Results.relSpkTms;
cStruct = pfs.Results.cStruct;
verbose = pfs.Results.verbose;
%% Validation of time windows
respWin = cStruct.Response_window_s;
sponWin = -flip(respWin);
tmLm = cStruct.Viewing_window_s;
binWdth = cStruct.BinSize_s;
respDelta = diff(respWin);
sponDelta = diff(sponWin);
tmDelta = diff(tmLm);
if respDelta > tmDelta || sponDelta > tmDelta && verbose
    fprintf(2, 'Total time limits (%.2f - %.2f ms) enclose ', tmLm*1e3)
    fprintf(2, 'a smaller period than the response time!\n')
end
if any(respWin > tmLm(2)) || any(respWin < tmLm(1)) && verbose
    fprintf(2, '%.2f - %.2f ms given as response window\n', respWin*1e3)
    fprintf(2, '%.2f - %.2f ms given as time limits...\n', tmLm*1e3)
end
if any(sponWin > tmLm(2)) || any(sponWin < tmLm(1)) && verbose
    fprintf(2, '%.2f - %.2f ms given as spontaneous window\n', sponWin*1e3)
    fprintf(2, '%.2f - %.2f ms given as time limits...\n', tmLm*1e3)
end
if respDelta ~= sponDelta && verbose
    fprintf(2, '%.2f ms enclosed by response window and\n', respWin*1e3)
    fprintf(2, '%.2f ms enclosed by spontaneous window...\n', sponWin*1e3)
end
%% Binarizing the spiking activity associated to the response
psthOpts = {'BinWidth', binWdth, 'BinLimits', tmLm};
fnOpts = {'UniformOutput', false};
erOpts = {"ErrorHandler",@badSizeSamp};
[Ncl, Na] = size(relSpkTms(1).SpikeTimes);
M = [1, -1;0, 1];
getProps = @(x) round(M*(x'.^[0;1])*Na);
% Get a PSTH per cluster per trial
psthPerCluster_perTrial =...
    arrayfun(@(y) cellfun(@(x) histcounts(x, psthOpts{:}),...
    y.SpikeTimes, fnOpts{:}), relSpkTms, fnOpts{:});
psthtx = tmLm(1)+(binWdth/2):binWdth:tmLm(2)-(binWdth/2);
% Organise the trials for a cluster in a matrix size #Trials x #Bins
p_f =...
    cellfun(@(y) arrayfun(@(x) cat(1, y{x,:}), (1:Ncl)', fnOpts{:}),...
    psthPerCluster_perTrial, fnOpts{:});
% Define a function that gets the mean spiking activity in a given window
getSpikingActivityIn = @(x, tw) mean(x(:, psthtx >= tw(1) & psthtx <= tw(2)),2);
pSpon = cellfun(@(y) cellfun(@(x)...
    mean(getSpikingActivityIn(x,sponWin) > 0), y), p_f, fnOpts{:});
pEvok = cellfun(@(y) cellfun(@(x)...
    mean(getSpikingActivityIn(x,respWin) > 0), y), p_f, fnOpts{:});
[h, pfs] = cellfun(@(x,y) arrayfun(@(u,v) ...
    fishertest(table(getProps(u), getProps(v))), x, y, fnOpts{:}), ...
    pSpon, pEvok, fnOpts{:});
h = gatherValues(h); pfs = gatherValues(pfs);
pwr = cellfun(@(x,y) arrayfun(@(u,v) sampsizepwr('p', u, v, [], Na), ...
    x, y, erOpts{:}), pEvok, pSpon, fnOpts{:});
pwr = cat(2, pwr{:});
pEvok = cat(2, pEvok{:});
pSpon = cat(2, pSpon{:});
p_f = cellfun(@(y) cellfun(@(x)...
    mean(getSpikingActivityIn(x,respWin)...
    > getSpikingActivityIn(x,sponWin)), y), p_f, fnOpts{:});
p_f = cat(2, p_f{:});
fishStats = struct('H',h,'P',pfs,'Power',pwr);
end

function out = badSizeSamp(S, varargin)
out = NaN;
end

function x_out = gatherValues(x_in)
x_out = cellfun(@(x) cell2mat(x), x_in, 'UniformOutput', false);
x_out = cat(2, x_out{:});
end