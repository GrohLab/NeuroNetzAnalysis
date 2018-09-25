function [featuresPerMiniCluster, meanWfStack, spTmcl] =...
    getClusterDetails(rootNames, clusterID)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pckNr = convertStringsToChars(...
    extractBetween(convertCharsToStrings(clusterID),3,"_"));
pckFile = strcat(rootNames,'_',pckNr,'.mat');
endCell = @(x) x{end};
if exist(pckFile,'file')
    load(pckFile,'spikes')
    clNr = str2double(endCell(strsplit(clusterID,'_')));
    if exist('spikes','var') 
        clIdx = spikes.assigns == clNr;
        clWf = spikes.waveforms(clIdx,:,:);
        clPCAu = spikes.info.pca.u(clIdx,:);
        clPCAs = spikes.info.pca.s;
        clPCAv = spikes.info.pca.v;
        miniCls = unique(spikes.info.kmeans.assigns(clIdx));
        Nch = numel(spikes.info.detect.stds);
        meanWfStack = zeros(length(miniCls),size(spikes.waveforms,2),Nch,'single');
        featuresPerMiniCluster = zeros(sum(clIdx),3,'single');
        MCidx = zeros(sum(clIdx),1,'uint8');
        spTmcl = zeros(sum(clIdx)-length(miniCls),1,'single');
        cmsum = 1; 
        spT = spikes.spiketimes;
        isi = [spT(1),diff(spT)];
        for cmcl = 1:length(miniCls)
            mclIdx = spikes.info.kmeans.assigns == miniCls(cmcl);
            meanWfStack(cmcl,:,:) = mean(spikes.waveforms(mclIdx,:,:),1);
            featuresPerMiniCluster(cmsum:cmsum+sum(mclIdx)-1,:) =...
                spikes.waveforms(mclIdx,:) * clPCAv(:,1:3);
            isimcl = isi(mclIdx);
            binEdges = exp(linspace(log(min(isi)),log(max(isi)),128));
            invFig = figure('Visible','off');hisi = histogram(isimcl,binEdges);
            visFig = figure;semilogx(binEdges,[0,hisi.BinCounts]);close(invFig)
            title(['ISI for mini-cluster ',num2str(miniCls(cmcl))])
            xlabel('\deltat [s]');ylabel('Counts')
            spTmcl(cmsum:cmsum+sum(mclIdx)-1) = isimcl;
            MCidx(cmsum:cmsum+sum(mclIdx)-1) = miniCls(cmcl);
            cmsum = cmsum + sum(mclIdx);
        end
        featuresPerMiniCluster = [featuresPerMiniCluster,MCidx];
    else
    end
    
else
end
for cf = 1:length(miniCls)
figure
for cch = 1:Nch
plot(meanWfStack(cf,:,cch))
if cch == 1
hold on
end
title(['Mini-cluster ',num2str(miniCls(cf))])
end
end
end

