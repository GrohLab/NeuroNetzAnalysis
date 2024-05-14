
fnOpts = {'UniformOutput', false};
[Ncl, Ntr] = cellfun(@size, {relativeSpkTmsStruct(:).SpikeTimes} );
Nccond = size( relativeSpkTmsStruct, 2 );
cwidth = 0.01; time_step = 0.001;
vWin = configStructure.Viewing_window_s;
bin_size = configStructure.BinSize_s;

[bC, bE] = prepareLogBinEdges([1e-5, cwidth], 63);
non_alloc_cells = cell( Nccond, max(Ntr) );
for cd = 1:Nccond
    fprintf(1, "Condition %s\n", relativeSpkTmsStruct(cd).name)
    for tr = 1:size(relativeSpkTmsStruct(cd).SpikeTimes, 2)
        fprintf(1, "Trial: %d\n", tr)
        rsp_ct = cat( 2, relativeSpkTmsStruct(cd).SpikeTimes{:,tr} );
        cw = vWin(1);
        while cw + time_step <= vWin(2)
            all_cl_spkTms = rsp_ct( rsp_ct > cw & ...
                rsp_ct < (cw + cwidth));
            if numel( all_cl_spkTms ) < 2
                intercl_isi = [];
            else    
                triu_subs = nchoosek( 1:numel( all_cl_spkTms ), 2 );
                intercl_isi = distmatrix( all_cl_spkTms(:), ...
                    all_cl_spkTms(:), 2 );
                intercl_isi = arrayfun(@(x,y) intercl_isi(x,y), ...
                    triu_subs(:,1), triu_subs(:,2) );
            end
            non_alloc_cells{cd,tr} = [non_alloc_cells{cd,tr}; 
                histcounts(intercl_isi , "BinEdges", [0,10.^bE] )];
            cw = cw + time_step;
        end
    end
end

time_isi_dist = arrayfun(@(c) cat( 3, non_alloc_cells{c,:} ), 1:Nccond, ...
    fnOpts{:} );
mdl_tx = fit_poly( [1, size(time_isi_dist{1}, 1)], ...
    vWin + [1,-1]*(time_step/2), 1 );
tx = ( (1:size(time_isi_dist{1}, 1) )'.^[1,0] ) * mdl_tx;
%%
figure_dir = ...
    "Z:\Jesus\Jittering\FULLYCurated\14_190716-2_jittering_3725_1620_1500 yes VPM good\KS2_newChanMap\Figures\Ephys VW-300.00-150.00 RW2.00-12.00 SW-212.00--202.00";
consCondNames = {relativeSpkTmsStruct(:).name};
for cc = 1:Nccond
    fig = figure(cc); imagesc(  tx,  bC , mean( time_isi_dist{cc}, 3)' ) 
    yticklabels( 1e3 * 10.^yticks ); xticklabels( 1e3 * xticks )
    ylabel('ISI [ms]'); xlabel('Time [ms]')
    title(['Intercluster ISI time distribution'; ...
        sprintf("%s", relativeSpkTmsStruct(cc).name) ] )
    clim([0, 19]); colormap(inferno)
    colorbar('Box', 'off')
    saveFigure( fig, fullfile( figure_dir, join([ ...
        "Intercluster ISI time dist", string(consCondNames(cc))]) ), true )
end

