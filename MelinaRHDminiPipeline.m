read_Intan_RHD2000_file
fs = frequency_parameters.amplifier_sample_rate;
filteredSignals = iir50NotchFilter(amplifier_data',fs);
filteredSignals = iirSpikeFilter(filteredSignals,fs)';
[Nch, Ns] = size(filteredSignals);
clearvars -except filteredSignals fs baseName path Nch Ns
%% 
spkSrt = UMSDataLoader(filteredSignals',fs); % Spike sort
spkSrt.UMS2kPipeline
spkStr = spkSrt.SpikeUMSStruct; % Spike structure from UMS2k
%% 
st = cell(1,sum(spkStr.labels(:,2) < 4));
ccl = 1;
ccl2 = ccl;
while ccl <= length(st)
    if spkStr.labels(ccl2,2) < 4
        st(ccl) = {spkStr.spiketimes(spkStr.assigns == spkStr.labels(ccl2,1))};
        ccl = ccl + 1;
    end
    ccl2 = ccl2 + 1;
end
%% Plot
labels = cell(Nch,1);
for cs = 1:Nch
    labels(cs) = {num2str(cs)};
end

%%
nch = 1:Nch;
%%
pltDur = 2;
sclFact = 5*numel(nch)*std(filteredSignals(:));
iSub = randi(Ns-fs,1);
fSub = iSub + (pltDur * 3e4);
plotEEGchannels(filteredSignals(nch,iSub:fSub),labels(nch),pltDur,fs,0.5);
fprintf('Showing from %.3f seconds\n',iSub/fs)

Npn = length(st);
hold on;
cmap = jet(numel(st));
for cpn = 1:Npn
    shwnSpk = iSub/fs <= st{cpn} & fSub/fs >= st{cpn};
    Nsk = sum(shwnSpk);
    plot(repmat(st{cpn}(shwnSpk),2,1)-(iSub/fs),...
        [sclFact*ones(1,Nsk,'single');zeros(1,Nsk,'single')],...
        'Color',cmap(cpn,:),'LineStyle','--','LineWidth',0.1)
    plot(st{cpn}(shwnSpk)-(iSub/fs),cpn*ones(1,Nsk,'single'),...
        'LineStyle','none','Marker','.','Color',cmap(cpn,:),...
        'DisplayName',sprintf('P. neuron %d',cpn));
end

%% Unit overlay

ccl2 = 1;
Nts = size(spkStr.waveforms,2);
tx = (0:Nts-1)*(1/fs);
f = gobjects(size(spkStr.labels,1),1);
for ccl = 1:size(spkStr.labels,1)
    if spkStr.labels(ccl,2) < 4
        spSub = spkStr.labels(ccl,1);
        spIdx = spkStr.assigns == spSub;
        swf_allChannels = spkStr.waveforms(spIdx,:,:);
        avSpkWvf = reshape(mean(swf_allChannels,1),Nts,Nch);
        [~, pch] = max(sum(abs(avSpkWvf)));
        pWvf = swf_allChannels(:,:,pch);
        f(ccl2) = figure('Name',sprintf('Cluster %d, channel %d',spSub,pch),...
            'Color',[1,1,1]);
        p = plot(tx,pWvf','LineWidth',0.1,'Color',[0.8,0.8,0.8]);box off;
        hold on;ap = plot(tx,avSpkWvf(:,pch),'LineWidth',2,'Color',[0,0,0]);
        legend([p(end),ap],{'Raw traces','Average trace'});
        title(sprintf('Cluster %d, channel %d',spSub,pch))
%         f(ccl2).PaperOrientation = 'landscape';
%         print(f(ccl2),fullfile(path,...
%             [baseName,sprintf('_cluster_%d.pdf',spSub)]),'-dpdf',...
%             '-bestfit')
%         close(f(ccl2));
        ccl2 = ccl2 + 1;
    end
end