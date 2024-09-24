function [mdlAll_ind, params] = regressEphysVSBehaviour( data_path, params )

fnOpts = {'UniformOutput', false};
tocol = @(x) x(:);
getAbsPath = @(x) string( fullfile( x.folder, x.name ) );

eph_path = dir( fullfile( data_path, "ephys*" ) );
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
for x=bfs_paths, load( getAbsPath( x ) ), end

eph_pttrns = ["*_Spike_Times.mat", "*analysis.mat"];
efs_paths = arrayfun(@(pt) dir( fullfile( eph_path, pt) ), eph_pttrns );
if any( ~arrayfun(@(x) exist( getAbsPath(x), "file" ), efs_paths ) )
    fprintf(1, 'Not all necessary ephys files exist!\n')
    return
end
for x=efs_paths, load( getAbsPath( x ) ), end

stop_time = length( Triggers.Whisker )/ fs;
bp_names = ["Stim-whisker mean", "Stim-whisker fan arc", ...
    "Nonstim-whisker mean", "Nonstim-whisker fan arc", ...
    "Interwhisker arc", "Symmetry", "Nose", "Roller speed"];

analysis_pttrn = "CW%.2f-%.2fms DW%.2f-%.2f BZ%.2f";
%%
behSignals = [behDLCSignals, vf];
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
ttx = (rel_win(1) + bin_size/2):bin_size:(rel_win(2) - bin_size/2);
%%
bin_edges = 0:bin_size:stop_time;
bin_centres = mean( [bin_edges(1:end-1); bin_edges(2:end)] );
Ntb = length( bin_centres );
hstOpts = {'Normalization', 'countdensity'};
binned_spikes = cellfun(@(s) histcounts( s, bin_edges, hstOpts{:}), ...
    spike_times, fnOpts{:} );
binned_spikes = cat( 1, binned_spikes{:} );

binned_beh = zeros( Ntb, Ns );
parfor b = 1:Ntb
    idx = my_xor( btx < bin_edges(b:b+1) );
    binned_beh(b,:) = mean( behSignals( idx , : ), 1 );
end

%% Design matrix for a set of trials (Control)
ctrl_sub = ismember( string( {Conditions.name} ), "Control Puff" );
time_limits = Conditions(ctrl_sub).Triggers(:,1)./fs + rel_win;
Nr = size( time_limits, 1 );
Nd = ceil( diff( del_win ) / bin_size );
auX = zeros( Nb*Nr, Nu, Nd );

cwin = arrayfun(@(x) linspace( time_limits(x,1) + (bin_size/2), ...
    time_limits(x,2) - (bin_size/2), Nb )', (1:Nr)', fnOpts{:} );
cwin = cat( 1, cwin{:} );

bin_ax = cwin + linspace( del_win(1)+(bin_size/2), ...
    del_win(2)-(bin_size/2), Nd );
% tr_ID = ceil( ( 1:(Nr*Nb) )' / Nb );
parfor r = 1:(Nr*Nb)
    tempC = my_cat( arrayfun( @(u) interp1( bin_centres, binned_spikes(u,:), ...
        bin_ax(r,:) ), 1:Nu, fnOpts{:} ), 1);
    auX( r, :, :) = tempC;
end

X = reshape( auX, Nb*Nr, Nu*Nd ); clearvars auX;
Xp = [ ones( Nb*Nr, 1), X];

%% Delay
delay_sub = cellfun(@(x) ~isempty(x), regexp( string( {Conditions.name} ), ...
    'Delay \d\.\d+\ss\s\+\sL' ) );

time_limits = Conditions(delay_sub).Triggers(:,1)./fs + rel_win;
Nr = size( time_limits, 1 );
Nd = ceil( diff( del_win ) / bin_size );
auX = zeros( Nb*Nr, Nu, Nd );

cwin = arrayfun(@(x) linspace( time_limits(x,1) + (bin_size/2), ...
    time_limits(x,2) - (bin_size/2), Nb )', (1:Nr)', fnOpts{:} );
cwin = cat( 1, cwin{:} );

bin_ax = cwin + linspace( del_win(1)+(bin_size/2), ...
    del_win(2)-(bin_size/2), Nd );
% tr_ID = ceil( ( 1:(Nr*Nb) )' / Nb );
parfor r = 1:(Nr*Nb)
    tempC = my_cat( arrayfun( @(u) interp1( bin_centres, binned_spikes(u,:), ...
        bin_ax(r,:) ), 1:Nu, fnOpts{:} ), 1);
    auX( r, :, :) = tempC;
end

X = reshape( auX, Nb*Nr, Nu*Nd ); clearvars auX;
Xl = [ ones( Nb*Nr, 1), X];
%% Multivariate regression response matrix
%X2 = [ ones( Nb*Nr, 1), X];
%lmObjs = cell( Ns, 1 );
y = zeros( Nb*Nr, Ns);
for r = 1:Nr
    idx = (r-1)*Nb + (1:Nb);
    aux = arrayfun(@(s) interp1( bin_centres, binned_beh(:,s), ...
        (1:Nb)'*bin_size + time_limits(r,1) ), (1:Ns), fnOpts{:} );
    aux = cat( 2, aux{:} );
    y(idx,:) = aux;
end
%{
%% Linear regression using fitlm
cvk = 15;
tr_ID = tocol( ones( Nb, 1 ) * (1:Nr) );
rmse1v = zeros( cvk , 1 ); mdl1v = cell( size( rmse1v ) );
Nk = round( Nr*0.15 ); 
parfor ii = 1:cvk
    testTrials = sort( randperm( Nr, Nk ) );
    trainingTrials = setdiff( 1:Nr, testTrials );
    trainingIdx = any( tr_ID == trainingTrials(:)', 2 );
    testIdx = ~trainingIdx;

    mdl1v{ii} = fitlm( X(trainingIdx,:), y(trainingIdx,1) );
    y_pred = predict( mdl1v{ii}, X(testIdx,:) );
    rmse1v(ii) = sqrt( mean( ( y(testIdx,1) - y_pred ).^2 ) );
end

[~, min_error] = min(rmse1v);
y_1_pred = predict( mdl1v{min_error}, X );

y_1 = reshape( y(:,1), Nb, Nr );
y_1_pred = reshape( y_1_pred, Nb, Nr );

%% Linear regression using ridge regularisation
cvk = 15; Nlambda = 64;
tr_ID = tocol( ones( Nb, 1 ) * (1:Nr) );
rmse1v = zeros( cvk , Nlambda ); 
mdl1v = zeros( size(X,2)+1, Nlambda, cvk );
Nk = round( Nr*0.15 ); 
lambdas = logspace( 3, 8, Nlambda);
parfor ii = 1:cvk
    testTrials = sort( randperm( Nr, Nk ) );
    trainingTrials = setdiff( 1:Nr, testTrials );
    trainingIdx = any( tr_ID == trainingTrials(:)', 2 );
    testIdx = ~trainingIdx;
    
    Xtrain = [ones( sum( trainingIdx ), 1 ), X(trainingIdx,:)];
    Xtest = [ones( sum( testIdx ), 1 ), X(testIdx,:)];
    % theta_0 = ( X2' * X2 ) \ ( X2' * y(trainingIdx,1) ) ;
    mdl1v(:,:,ii) = ridge( y(trainingIdx,1), Xtrain, lambdas );
    y_pred = Xtest * mdl1v(:,:,ii);
    rmse1v(ii,:) = sqrt( mean( ( y(testIdx,1) - y_pred ).^2 ) );
end

%}
%% Linear regression for all behavioural signals using matrix multiplication
cvk = 15;
tr_ID = tocol( ones( Nb, 1 ) * (1:Nr) );
rmseAll_ind = zeros( cvk, Ns );
mdlAll_ind = zeros( size(X,2)+1, cvk, Ns );
Nk = round( Nr*0.15 ); idxs = zeros( cvk, Nk );
parfor ii = 1:cvk
    testTrials = sort( randperm( Nr, Nk ) );
    idxs(ii,:) = testTrials;
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
analysis_key = sprintf( analysis_pttrn, rel_win/m, del_win/m, bin_size/m );
save( fullfile( data_path,  join( ["Regression", analysis_key + ".mat"] ) ), ...
    "-v7.3" )
%%
createtiles = @(f,nr,nc) tiledlayout( f, nr, nc, ...
    'TileSpacing', 'Compact', 'Padding', 'tight');
cleanAxis = @(x) set( x, "Box", "off", "Color", "none" );
figWeight = figure('Color', 'w'); t = createtiles( figWeight, 2, 4 );
mdlAll_ind_norm = [mdlAll_ind(1,:,:);
    mdlAll_ind(2:end,:,:) ./ vecnorm( mdlAll_ind(2:end,:,:), 2, 1 )];
mdl_mu = squeeze( mean( mdlAll_ind_norm, 2 ) );
for cb = 1:Ns
    ax = nexttile(t);
    imagesc( ax, wtx/m, [], reshape( mdl_mu(2:end,cb), Nu, Nd ) )
    cleanAxis( ax ); yticks( ax, 1:Nu ); title( ax, bp_names( cb ) );
    colormap( traffic ); 
    clim( 1.3*max( abs( mdl_mu(2:end,cb) ) )*[-1,1] )
    cbObj = colorbar( 'Box', 'off', 'AxisLocation', 'out', ...
        'TickDirection', 'out', 'Location', 'northoutside' );
end
xlabel(ax, 'Time [ms]'); axs = findobj( t, "Type", "Axes" );
ylabel( axs(end), 'Units' )
title( t, 'Regression weights' )
arrayfun(@(x) set( get( x, "YAxis" ), "Visible", "off" ), ...
    axs(setdiff( 1:Ns, [4,8] )) )
arrayfun(@(x) set( get( x, "XAxis" ), "Visible", "off" ), axs(5:8) )

saveFigure( figWeight, fullfile( eph_path, "Figures", ...
    join( ["Regression weights", analysis_key] ) ), true )
%%
errFig = figure("color", "w"); t = createtiles( errFig, 1, 1 );
ax = nexttile(t);
gray15pc = 0.15*ones(1,3);
boxchart(ax, rmseAll_ind./range(y) , 'Notch', 'on', ...
    'BoxFaceColor', gray15pc, 'JitterOutliers', 'on', ...
    'MarkerStyle', '.', 'MarkerColor', gray15pc );
xticklabels( ax, bp_names ); cleanAxis( ax ); 
ylabel( ax, 'Normalised error' )
title( ax, sprintf( '%d-fold cross-validated error', cvk ) )

saveFigure( errFig, fullfile( eph_path, "Figures", ...
    join( [sprintf( "%d-fold cv error", cvk ), analysis_key] ) ), true )
%{
%% Multiple output linear regression
cvk = 15;
Nk = round( Nr*0.15 );
rmse = zeros( cvk , Ns ); 
mdl = zeros( size( X, 2 )+Ns, size( y, 2 ), cvk );
idxs = zeros( cvk, Nk ); %[zy, y_mu, y_sig] = zscore(y, 0, 1);
X2 = [[eye(Ns); zeros( size(X,1) - Ns, Ns )], X];
for ii = 1:cvk
    fprintf(1, 'K:%d\n', ii)
    testTrials = sort( randperm( Nr, Nk ) );
    idxs(ii,:) = testTrials;
    trainingTrials = setdiff( 1:Nr, testTrials );
    trainingIdx = any( tr_ID == trainingTrials, 2 );
    testIdx = ~trainingIdx;

    mdl(:,:,ii) = mvregress( gpuArray( X2(trainingIdx,:) ), ...
        gpuArray( y(trainingIdx,:) ) );
    y_pred = X2(testIdx,:) * mdl(:,:,ii);
    rmse(ii,:) = sqrt( mean( ( y(testIdx,:) - y_pred ).^2 ) );
end

[~, min_error] = min(rmse,[],1);
y_all_pred = X2 * squeeze( mean( mdl, 3 ) );

y_trials = reshape( y, Nb, Nr, Ns );
y_all_pred = reshape( y_all_pred, Nb, Nr, Ns );

%% Training

ho_trials = randperm( Nr, round( Nr*0.1 ) );
testIdx = any( tr_ID == sort(ho_trials), 2 );
cv_kf = cvpartition( tr_ID( ~testIdx ), "KFold", 15 );

Xtrain = X(~testIdx, :); ytrain = y(~testIdx,1);

[mdl, fitInfo] = lassoglm( X(~testIdx,:), y(~testIdx,1), 'normal', ...
    'CV', cv_kf, 'Lambda', logspace( -5, 3, 64 ), ...
    'Options', statset('UseParallel', true ), ...
    'Alpha', eps );

w_vec = [fitInfo.Intercept(fitInfo.IndexMinDeviance);
    mdl(:,fitInfo.IndexMinDeviance)];
y_pred = glmval( w_vec, X, "identity" );
y_pred = reshape( y_pred, Nb, Nr );
clrMap = [0.15*ones(1,3); 0.85,0.51,0.15 ];
figure; lObj = line( 1:(Nb*numel( ho_trials )), [y(testIdx,1), y_pred(:)] );
arrayfun(@(ii,x) set( x, 'Color', clrMap(ii,:) ), (1:numel(lObj))', lObj(:) )

createtiles = @(f,nr,nc) tiledlayout( f, nr, nc, ...
    'TileSpacing', 'Compact', 'Padding', 'tight');
cleanAxis = @(x) set( x, "Box", "off", "Color", "none" );

fig = figure( 'Color', 'w' ); t = createtiles( fig, 10, 2);
nexttile([8, 1]);
imagesc( ttx, [], y_1' - 90); colormap(inferno)
xline(0, 'LineStyle', '--', 'Color', 0.85*ones(1,3) )
nexttile([8, 1]);
imagesc( ttx, [], y_pred' - 90); colormap(inferno)
set( get( gca, 'YAxis' ), 'Visible', 'off' )
nexttile(t);
line( ttx, mean( y_1, 2 ) - 90, 'Color', 0.15*ones(1,3), 'LineWidth', 1.5 )
nexttile(t);
line( ttx, mean( y_pred, 2 ) - 90, 'Color', [0.85, 0.51, 0.15] , 'LineWidth', 1.5 )
axs = get( t, "Children" );
linkaxes( axs, 'x')
xlim(ttx([1,end]))
arrayfun(@(x) set( get( x, "XAxis" ), "Visible", "off" ), axs(3:4) )
arrayfun(@(x) xticklabels( x, xticks(x) / m ), axs(1:2) )

rmse = mean( ( y_1 - y_pred ).^2, 1 );


wtx = (del_win(1) + bin_size/2):bin_size:(del_win(2) - bin_size/2);
ttx = (rel_win(1) + bin_size/2):bin_size:(rel_win(2) - bin_size/2);

figure; imagesc( wtx, [], reshape( w_vec(2:end), Nd, Nu )' )

%% Design matrix for the whole experiment
Nb = ceil( diff( bin_edges([1,end]) )/ bin_size );
auX = zeros( Nb, Nu, Nd );
parfor b = 1:Nb
    cwin = bin_centres(b) + del_win;
    bin_ax = linspace( cwin(1), cwin(2), Nd );
    tempC = arrayfun(@(u) interp1( bin_centres, binned_spikes(u,:), ...
        bin_ax ), 1:Nu, fnOpts{:} );
    tempC = cat( 1, tempC{:} );
    tempC( isnan(tempC) ) = 0;
    auX( b, :, :) = tempC;
end
X = reshape( auX, [], Nu*Nd );
X2 = [ ones( Nb*Nr, 1), X];
Xa = X2;
%}
end