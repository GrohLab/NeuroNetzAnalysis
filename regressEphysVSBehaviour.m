function [mdlAll_ind, params, DX] = regressEphysVSBehaviour( data_path, params, varargin )
%% Input parsing
p = inputParser;
my_xor = @(x) xor( x(:,1), x(:,2) );

defCond = 'freq';
cont_exp = '^Delay\s\d\.\d+\ss';
freq_exp = [cont_exp, '\s\+\sL'];
checkCondSel = @(x) (isstring(x) | ischar(x)) & ...
    my_xor(string(x) == ["freq", "cont"]);
param_fields = {'relative_window', 'delay_window', 'bin_size', 'kfold'};

checkParams = @(x) isstruct( x ) & all( isfield( x, param_fields ) );

addParameter( p, 'Condition', defCond, checkCondSel )
addParameter( p, 'Verbose', true, @(x) islogical(x) & isscalar(x) )
addRequired( p, 'data_path', @(x) exist(data_path, "dir") )
addRequired( p, 'params', checkParams )

parse(p, data_path, params, varargin{:} )

selCond = p.Results.Condition;
verbose = p.Results.Verbose;
data_path = p.Results.data_path;
params = p.Results.params;

cond_exp = freq_exp;
if ~strcmp( selCond, "freq" )
    cond_exp = cont_exp;
end
%% 

fnOpts = {'UniformOutput', false};
tocol = @(x) x(:);
getAbsPath = @(x) string( fullfile( x.folder, x.name ) );
ovwFlag = false;
eph_path = dir( fullfile( data_path, "ephys*" ) );
mdlAll_ind = nan; DX = cell(4,1);
if ~isempty( eph_path )
    eph_path = getAbsPath( eph_path );
else
    fprintf(1, 'No ephys folder found!\n')
    return
end
beh_path = fullfile( data_path, "Behaviour" );

beh_pttrns = ["RollerSpeed*.mat", "BehaviourSignals*.mat"];
bfs_paths = arrayfun(@(pt) dir( fullfile( beh_path, pt) ), beh_pttrns );
if any( ~arrayfun(@(x) exist( getAbsPath(x), "file" ), bfs_paths ) )
    fprintf(1, 'Not all necessary behaviour files exist!\n')
    return
end
for x = bfs_paths, load( getAbsPath( x ) ), end

eph_pttrns = ["*_Spike_Times.mat", "*analysis.mat"];
try
    efs_paths = arrayfun(@(pt) dir( fullfile( eph_path, pt) ), eph_pttrns );
catch
    fprintf(1, "Seems like no in-depth analysis was run in this directory\n")
    fprintf(1, "%s\n", data_path )
    return
end
if any( ~arrayfun(@(x) exist( getAbsPath(x), "file" ), efs_paths ) )
    fprintf(1, 'Not all necessary ephys files exist!\n')
    return
end
for x = efs_paths, load( getAbsPath( x ) ), end

stop_time = length( Triggers.Whisker )/ fs;
bp_names = ["Stim-whisker mean", "Stim-whisker fan arc", ...
    "Nonstim-whisker mean", "Nonstim-whisker fan arc", ...
    "Interwhisker arc", "Symmetry", "Nose", "Roller speed"];

analysis_pttrn = "CW%.2f-%.2fms DW%.2f-%.2f BZ%.2f";

figOpts = {'Visible','on'};
ofgOpts = {'new', 'visible'};
if ~strcmp(computer, 'PCWIN64')
    figOpts{2} = 'off';
    ofgOpts{2} = 'invisible';
end
vars2save = {'mdlAll_ind', 'params', 'DX', 'analysis_key', ...
    'binned_spikes', 'binned_beh'};

delay_sub = cellfun(@(x) ~isempty(x), regexp( string( {Conditions.name} ), ...
    cond_exp ) );

if all(~delay_sub)
    fprintf( 1, 'No laser in frequency used!\nAborting.\n')
    return
end
%%
try
    behSignals = [behDLCSignals, vf];
catch
    fprintf(1, 'Dimensions of roller speed and other behaviour signals')
    fprintf(1, ' mismatch!\n')
    fprintf(1, "%s\n", data_path )
    return
end
mdl_btx = fit_poly( [1, size( behSignals, 1 )], [0, size( behSignals, 1 )/fr] + [1,-1] * (1/fr), 1 );
btx = (1:size( behSignals, 1 ))'.^[1,0] * mdl_btx;
my_xor = @(x) xor( x(:,1), x(:,2) );
my_cat = @(x,d) cat( d, x{:} );

m = 1e-3;
rel_win = params.relative_window; % [-1, 1]*0.8;
del_win = params.delay_window; % [-100, 100]*m;
bin_size = params.bin_size; % 10*m;

Nb = ceil( diff( rel_win )/ bin_size );
% Nb = ceil( diff( time_limits ) / bin_size );
Nu = numel( spike_times );
% cons_time = my_xor( btx > time_limits );
Ns = size( behSignals, 2 );
wtx = (del_win(1) + bin_size/2):bin_size:(del_win(2) - bin_size/2);
trial_tx = (rel_win(1) + bin_size/2):bin_size:(rel_win(2) - bin_size/2);

analysis_key = sprintf( analysis_pttrn, rel_win/m, del_win/m, bin_size/m );
regFile = fullfile( data_path,  join( ["Regression", analysis_key + ".mat"] ) );
rfeFlag = false;
if exist(regFile, "file")
    rfeFlag = true;
    load( regFile )
end
%%
bin_edges = 0:bin_size:stop_time;
bin_centres = mean( [bin_edges(1:end-1); bin_edges(2:end)] );
if ~exist( 'binned_beh', 'var' )
    Ntb = length( bin_centres );
    hstOpts = {'Normalization', 'countdensity'};
    binned_spikes = cellfun(@(s) histcounts( s, bin_edges, hstOpts{:}), ...
        spike_times, fnOpts{:} );
    binned_spikes = cat( 1, binned_spikes{:} );
    if verbose
        fprintf(1, "Binning behaviour signals\n")
    end
    binned_beh = zeros( Ntb, Ns );
    parfor b = 1:Ntb
        idx = my_xor( btx < bin_edges(b:b+1) );
        binned_beh(b,:) = mean( behSignals( idx , : ), 1 );
    end
end

%% Design matrix for a set of trials (Control)
ctrl_sub = ismember( string( {Conditions.name} ), "Control Puff" );
time_limits = Conditions(ctrl_sub).Triggers(:,1)./fs + rel_win;
Nr = size( time_limits, 1 );
Nd = ceil( diff( del_win ) / bin_size );
if ~exist( 'DX', 'var' ) || any( size( DX{2} ) ~= [Nb*Nr, 1 + Nu*Nd] )
    auX = zeros( Nb*Nr, Nu, Nd );

    cwin = my_cat( arrayfun(@(x) linspace( time_limits(x,1) + (bin_size/2), ...
        time_limits(x,2) - (bin_size/2), Nb )', (1:Nr)', fnOpts{:} ), 1 );

    bin_ax = cwin + linspace( del_win(1)+(bin_size/2), ...
        del_win(2)-(bin_size/2), Nd );
    if verbose
        fprintf(1, "Design matrix\n")
    end

    parfor r = 1:(Nr*Nb)
        tempC = my_cat( arrayfun( @(u) interp1( bin_centres, binned_spikes(u,:), ...
            bin_ax(r,:) ), 1:Nu, 'UniformOutput', false ), 1);
        auX( r, :, :) = tempC;
    end

    X = reshape( auX, Nb*Nr, Nu*Nd ); clearvars auX;
    Xp = [ ones( Nb*Nr, 1), X];
    DX{2} = Xp;
else
    Xp = DX{2};
end

params.Nb = Nb; params.Nr = Nr; params.Ns = Ns; params.Nd = Nd;
params.Nu = Nu;
%% Multivariate regression response matrix
if ~exist( 'DX', 'var' ) || any( size( DX{1} ) ~= [Nb*Nr, Ns] )
    y = zeros( Nb*Nr, Ns);
    for r = 1:Nr
        idx = (r-1)*Nb + (1:Nb);
        aux = arrayfun(@(s) interp1( bin_centres, binned_beh(:,s), ...
            (1:Nb)'*bin_size + time_limits(r,1) ), (1:Ns), fnOpts{:} );
        aux = cat( 2, aux{:} );
        y(idx,:) = aux;
    end
    DX{1} = y;
else
    y = DX{1};
end

%% Linear regression for all behavioural signals using matrix multiplication
if ~exist('mdlAll_ind', 'var') || ...
        numel(mdlAll_ind)<prod([1 + Nu*Nd, params.kfold, Ns]) %|| ...
        %any( size( mdlAll_ind ) ~= [1 + Nu*Nd, params.kfold, Ns] )
    cvk = params.kfold;
    if Nr <= cvk
        cvk = round( Nr*0.75 );
        params.kfold = cvk;
    end
    cvObj = cvpartition( Nr, "KFold", cvk );
    tr_ID = tocol( ones( Nb, 1 ) * (1:Nr) );
    rmseAll_ind = zeros( cvk, Ns );
    mdlAll_ind = zeros( size(Xp,2), cvk, Ns );
    % Nk = round( Nr*0.15 );
    if verbose
        fprintf(1, "Regression\n")
    end
    parfor ii = 1:cvk
        if verbose
            fprintf(1, "%d ", ii)
        end
        testTrials = find( test( cvObj, ii ) );
        trainingTrials = setdiff( 1:Nr, testTrials );
        trainingIdx = any( tr_ID == trainingTrials, 2 );
        testIdx = ~trainingIdx;

        Xtrain = Xp(trainingIdx,:); Xtest = Xp(testIdx,:);
        for cb = 1:Ns
            ytrain = y(trainingIdx, cb); ytest = y(testIdx, cb);
            mdlAll_ind(:,ii,cb) = ( Xtrain' * Xtrain ) \ ( Xtrain' * ytrain ) ;
            %mdl1v(:,:,ii) = ridge( y(trainingIdx,1), Xtrain, lambdas );
            y_pred = Xtest * mdlAll_ind(:,ii,cb);
            rmseAll_ind(ii,cb) = sqrt( mean( ( ytest - y_pred ).^2 ) );
        end
    end
    if verbose
        fprintf(1, "\n")
    end
    if (sum(isnan(mdlAll_ind),"all") / numel( mdlAll_ind ))  > 0.1
        fprintf(1, 'This regression is not possible. Many NaN values returned\n')
        return
    end
    params.fit_error = rmseAll_ind;
else
    rmseAll_ind = params.fit_error;
    cvk = params.kfold;
end
%%
if verbose
    fprintf(1, "Plotting\n")
end
createtiles = @(f,nr,nc) tiledlayout( f, nr, nc, ...
    'TileSpacing', 'Compact', 'Padding', 'tight');
cleanAxis = @(x) set( x, "Box", "off", "Color", "none" );
fwPath = fullfile( eph_path, "Figures", ...
        join( ["Regression weights", analysis_key] ) );
if ~exist( fwPath + ".fig" , "file" )
    figWeight = figure('Color', 'w', figOpts{:});
    t = createtiles( figWeight, 2, 4 );
    mdlAll_ind_norm = [mdlAll_ind(1,:,:);
        mdlAll_ind(2:end,:,:) ./ vecnorm( mdlAll_ind(2:end,:,:), 2, 1 )];
    mdl_mu = squeeze( mean( mdlAll_ind_norm, 2 ) );
    for cb = 1:Ns
        ax = nexttile(t);
        imagesc( ax, wtx/m, [], reshape( mdl_mu(2:end, cb), Nu, Nd ) )
        cleanAxis( ax ); yticks( ax, 1:Nu ); title( ax, bp_names( cb ) );
        colormap( traffic );
        clim( 1.3*max( abs( mdl_mu(2:end,cb) ) )*[-1,1] )
        colorbar( 'Box', 'off', 'AxisLocation', 'out', ...
            'TickDirection', 'out', 'Location', 'northoutside' );
    end
    xlabel(ax, 'Time [ms]'); axs = findobj( t, "Type", "Axes" );
    ylabel( axs(end), 'Units' )
    title( t, 'Regression weights' )
    arrayfun(@(x) set( get( x, "YAxis" ), "Visible", "off" ), ...
        axs(setdiff( 1:Ns, [4,8] )) )
    arrayfun(@(x) set( get( x, "XAxis" ), "Visible", "off" ), axs(5:8) )

    saveFigure( figWeight, fwPath, true, ovwFlag )
else
    openfig( fwPath + ".fig", ofgOpts{:} );
end
%%
efPath = fullfile( eph_path, "Figures", ...
    join( [sprintf( "%d-fold cv error", cvk ), analysis_key] ) );
if ~exist( efPath + ".fig", 'file' )
    errFig = figure("color", "w"); t = createtiles( errFig, 1, 1 );
    ax = nexttile(t);
    gray15pc = 0.15*ones(1,3);
    boxchart(ax, rmseAll_ind./range(y) , 'Notch', 'on', ...
        'BoxFaceColor', gray15pc, 'JitterOutliers', 'on', ...
        'MarkerStyle', '.', 'MarkerColor', gray15pc );
    xticklabels( ax, bp_names ); cleanAxis( ax );
    ylabel( ax, 'Normalised error' )
    title( ax, sprintf( '%d-fold cross-validated error', cvk ) )

    saveFigure( errFig, efPath, true, ovwFlag )
else
    openfig( efPath + ".fig", ofgOpts{:} );
end
%% Delay

time_limits = Conditions(delay_sub).Triggers(:,1)./fs + rel_win;
Nr = size( time_limits, 1 );
Nd = ceil( diff( del_win ) / bin_size );
if ~exist( 'DX', 'var' ) || any( size( DX{3} ) ~= [Nb*Nr, 1 + Nu*Nd] )

    auX = zeros( Nb*Nr, Nu, Nd );

    cwin = arrayfun(@(x) linspace( time_limits(x,1) + (bin_size/2), ...
        time_limits(x,2) - (bin_size/2), Nb )', (1:Nr)', fnOpts{:} );
    cwin = cat( 1, cwin{:} );

    bin_ax = cwin + linspace( del_win(1)+(bin_size/2), ...
        del_win(2)-(bin_size/2), Nd );
    if verbose
        fprintf(1, "Design matrix 2\n")
    end
    % tr_ID = ceil( ( 1:(Nr*Nb) )' / Nb );
    parfor r = 1:(Nr*Nb)
        tempC = my_cat( arrayfun( @(u) interp1( bin_centres, binned_spikes(u,:), ...
            bin_ax(r,:) ), 1:Nu, fnOpts{:} ), 1);
        auX( r, :, :) = tempC;
    end

    X = reshape( auX, Nb*Nr, Nu*Nd ); clearvars auX;
    Xl = [ ones( Nb*Nr, 1), X];
    DX{3} = Xl;
    params.Laser = struct( 'Nr', Nr, 'Nd', Nd );
% else
%     Xl = DX{3};
end

%% Multivariate regression response matrix (Laser)
ylFlag = false;
if ~exist( 'DX', 'var' ) || numel(DX)<4 || any( size( DX{4} ) ~= [Nb*Nr, Ns] )

    yl = zeros( Nb*Nr, Ns);
    for r = 1:Nr
        idx = (r-1)*Nb + (1:Nb);
        aux = arrayfun(@(s) interp1( bin_centres, binned_beh(:,s), ...
            (1:Nb)'*bin_size + time_limits(r,1) ), (1:Ns), fnOpts{:} );
        aux = cat( 2, aux{:} );
        yl(idx,:) = aux;
    end
    DX{4} = yl;
% else
%     yl = DX{4};
%     ylFlag = true;
end
%%
clearvars *fig* ax* cbObj 
% DX = {y, Xp, Xl, yl};
if ~rfeFlag
    save( regFile, vars2save{:}, "-v7.3" )
else
    vrs = whos('-file', regFile);
    vrsInFile = string({vrs(:).name});
    if any( ~contains( vrsInFile, vars2save ) ) || ...
            prod(vrs(arrayfun(@(x) strcmp( 'DX', x.name ), vrs )).size) < numel(DX)
        save( regFile, vars2save{:}, '-append' )
    end
end
end