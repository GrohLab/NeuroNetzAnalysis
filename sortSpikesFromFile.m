function [spikes,outputArg2] =...
    sortSpikesFromFile(spikes,inputFile,chanOrder,chPp,ovrlap)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Loading data
if exist(inputFile,'file')
    fid = fopen(inputFile);
    [fDir,baseName,fExt] = fileparts(fopen(fid));
    if strcmp(fExt,'.smrx') || strcmp(fExt,'.smr')
        currDir = cd;
        importSMR([baseName,fExt],fDir,-1)
        cd(currDir)
    end
    if ~exist('chanOrder','var') || isempty(chanOrder)
        chanOrder = [8, 9, 7, 10, 4, 13, 5, 12, 2, 15, 1, 16, 6, 11, 3, 14];
    end
    dataLoader = UMSDataLoader(fullfile(fDir,[baseName,'.mat']),chanOrder);
    %% Sorting the spikes
    data = dataLoader.getDataMatrix;
    Nch = dataLoader.Nch;
    % Channels per pack
    cch = 1;
    spikes.params.Fs = dataLoader.SamplingFrequency;
    pckCounter = 1;
    while cch + chPp - 1 <= Nch
        auxData{1,1} = data(:,cch:cch + chPp - 1);
%         spikes = ss_detect(auxData,spikes);
%         spikes = ss_align(spikes);
%         spikes = ss_kmeans(spikes);
%         spikes = ss_energy(spikes);
%         spikes = ss_aggregate(spikes);
%         splitmerge_tool(spikes)
%         spikes = waitForUserAndGetSpikesStructure();
%         outlier_tool(spikes)
%         spikes = waitForUserAndGetSpikesStructure();
%         save(fullfile(fDir,...
%             [baseName,'_Pack',num2str(pckCounter),'.mat']),'spikes')
        pckCounter = pckCounter + 1;
        cch = cch + chPp - ovrlap;
    end
end
outputArg2 = chanOrder;
end

