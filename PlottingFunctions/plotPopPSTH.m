function ax = plotPopPSTH(PSTHStruct, varargin)
%PLOTPOPPSTH This function accepts the PSTH structure created by getPopPSTH
%and plots all PSTHs of all experiments together with a population PSTH at
%the bottom of the figure. It can plot the PSTH structure in a given axis
%or create a new figure to plot this results.
%   Emilio Isaías-Camacho @ GrohLab 2019

%% Input validation
if isstruct(PSTHStruct) && isfield(PSTHStruct,'TimeAxis')
    p = inputParser;
    defaultFigure = 'new';
    validFigure = {'new','current','given'};
    checkFigure = @(x) any(validatestring(x,validFigure));
    
    addRequired(p,'PSTHStruct',@isstruct);
    defaultBin = PSTHStruct.BinSize;
    small_limit = mean(diff(PSTHStruct.TimeAxis));
    attributes = {'scalar','>',small_limit,...
        '<=', sum(abs([PSTHStruct.TimeAxis(1),PSTHStruct.TimeAxis(end)]))*0.01};
    addOptional(p,'BinSize', defaultBin, @(x)validateattributes(x,...
        {'numeric'},attributes));
    addOptional(p,'Figure',defaultFigure, checkFigure);
    p.parse(PSTHStruct, varargin{:})
else
    fprintf('The PSTH structure is either corrupt or hasn''t been created\n')
    return
end
%% Initialization
binsize = p.Results.BinSize;
% Number of conditions combinations
Ncc = length(PSTHStruct.PSTHs);
% Figure organization
if sqrt(Ncc) - fix(sqrt(Ncc)) == 0
    Nrow = sqrt(Ncc);
    Ncol = Nrow;
elseif ~mod(Ncc,2)
    Nrow = 2;
    Ncol = Ncc/Nrow;
elseif ~mod(Ncc,3)
    Nrow = 3;
    Ncol = Ncc/Nrow;
end


end

% p = inputParser;
% charchk = {'char'};
% numchk = {'numeric'};
% nempty = {'nonempty'};
%
% addRequired(p,'shape',@(x)validateattributes(x,charchk,nempty))
% addRequired(p,'dim1',@(x)validateattributes(x,numchk,nempty))
% addOptional(p,'dim2',1,@(x)validateattributes(x,numchk,nempty))
% parse(p,shape,dim1,varargin{:})
