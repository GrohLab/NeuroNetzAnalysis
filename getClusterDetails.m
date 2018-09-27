function [featuresPerMiniCluster, miniIdx, ISImcl] =...
    getClusterDetails(rootNames, clusterID, plotFlag)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
featuresPerMiniCluster = 0;
miniIdx = false;
ISImcl = 0;
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = true;
end

pckNr = convertStringsToChars(...
    extractBetween(convertCharsToStrings(clusterID),3,"_"));
pckFile = strcat(rootNames,'_',pckNr,'.mat');
endCell = @(x) x{end};
if exist(pckFile,'file')
    load(pckFile,'spikes')
    clNr = str2double(endCell(strsplit(clusterID,'_')));
    if exist('spikes','var')
        clIdx = spikes.assigns == clNr;
        clPCAv = spikes.info.pca.v;
        miniCls = unique(spikes.info.kmeans.assigns(clIdx));
        Nch = numel(spikes.info.detect.stds);
        meanWfStack = zeros(length(miniCls),size(spikes.waveforms,2),Nch,'single');
        miniIdx = false(length(miniCls),length(clIdx));
        featuresPerMiniCluster = zeros(sum(clIdx),3,'single');
        MCidx = zeros(sum(clIdx),1,'uint8');
        ISImcl = zeros(sum(clIdx)-length(miniCls),1,'single');
        cmsum = 1;
        spT = spikes.spiketimes;
        isi = [spT(1),diff(spT)];
        tx = 0:1/spikes.params.Fs:(size(spikes.waveforms,2) - 1)/spikes.params.Fs;
        binEdges = exp(linspace(log(min(isi)),log(max(isi)),512));
        Ncl = numel(miniCls);
        cmap = jet(Ncl);
        figTimeVsISI = figure('Name','Time VS ISI','Color',[1,1,1]);
        for cmcl = 1:Ncl
            mclIdx = spikes.info.kmeans.assigns == miniCls(cmcl);
            miniIdx(cmcl,:) = mclIdx;
            meanWfStack(cmcl,:,:) = mean(spikes.waveforms(mclIdx,:,:),1);
            featuresPerMiniCluster(cmsum:cmsum+sum(mclIdx)-1,:) =...
                spikes.waveforms(mclIdx,:) * clPCAv(:,1:3);
            isimcl = isi(mclIdx);
            invFig = figure('Visible','off');hisi = histogram(isimcl,binEdges);
            figMC = figure('Name',['Mini-cluster ',num2str(miniCls(cmcl))],...
                'Color',[1,1,1]);
            axMC = subplot(2,4,[1,2,5,6]);
            figure(figTimeVsISI);
            axISI = subplot(1,4,1);semilogy(axISI,[0,hisi.BinCounts],...
                binEdges*1e3,'Color',cmap(cmcl,:),...
                'DisplayName',['ISI ',miniCls(cmcl)]);
            axtvi = subplot(1,4,2:4);semilogy(axtvi,...
                spT(miniIdx(cmcl,:)),isi(miniIdx(cmcl,:))*1e3,...
                'Color',cmap(cmcl,:),'DisplayName',num2str(miniCls(cmcl)),...
                'Marker','.','LineStyle','none');
            if cmcl == 1
                hold(axISI,'on')
                hold(axtvi,'on')
            end
            semilogx(axMC,binEdges*1e3,([0,hisi.BinCounts]),'Color',cmap(cmcl,:));
            close(invFig);title(axMC,['ISI for mini-cluster ',num2str(miniCls(cmcl))])
            xlabel(axMC,'$\frac{\delta}{\delta t} \vec{T}$ [ms]',...
                'Interpreter','latex','FontSize',13);ylabel(axMC,'Counts')
            subPlAx = [3,4,7,8];
            figure(figMC)
            for cch = 1:Nch
                ax = subplot(2,4,subPlAx(cch));
                plot(ax,tx,squeeze(meanWfStack(cmcl,:,cch)),'DisplayName',...
                    ['Channel ',num2str(cmcl)],'Color',cmap(cmcl,:));
            end
            
            ISImcl(cmsum:cmsum+sum(mclIdx)-1) = isimcl;
            MCidx(cmsum:cmsum+sum(mclIdx)-1) = miniCls(cmcl);
            cmsum = cmsum + sum(mclIdx);
        end
        set(axISI,'XDir','reverse','YGrid','on');xlabel(axISI,'Counts');
        ylabel(axISI,'$\frac{\delta}{\delta t} \vec{T}$ [ms]',...
                'Interpreter','latex','FontSize',13)
        set(axtvi,'YTickLabel',[],'YGrid','on');xlabel(axtvi,'Time [s]')
        title(axtvi,[clusterID,' Time VS ISI distribution'])
        linkaxes([axtvi,axISI],'y');legend(axtvi,'show','Location','best')
    else
        disp('The UMS ''spikes'' structure couldn''t be found in the corresponding file.')
        return
    end
else
    disp('Verify the rootName. Pack files couldn''t be found.')
    return
end
%% PLOTTING Features
figure('Name',['Features for ',clusterID],'Color',[1,1,1]);
for cfd = 1:Ncl
    scatter3(featuresPerMiniCluster(MCidx == miniCls(cfd),1),...
        featuresPerMiniCluster(MCidx == miniCls(cfd),2),...
        featuresPerMiniCluster(MCidx == miniCls(cfd),3),...
        15,cmap(cfd,:),'DisplayName',num2str(miniCls(cfd)),...
        'MarkerFaceColor',cmap(cfd,:),'MarkerEdgeAlpha',0.4,...
        'MarkerFaceAlpha',0.4)
    if cfd == 1
        hold on
    end
end
legend show;xlabel('PC 1');ylabel('PC 2');zlabel('PC 3');
title('Feature space for the selected clusters')
allFeat = spikes.waveforms(:,:) * clPCAv(:,1:3);
scatter3(allFeat(~clIdx,1),allFeat(~clIdx,2),allFeat(~clIdx,3),...
        5,[0.9,0.9,0.9],'DisplayName','Other units','MarkerEdgeAlpha',0.7,...
        'MarkerFaceColor',[0.9,0.9,0.9],'MarkerFaceAlpha',0.7)
end

