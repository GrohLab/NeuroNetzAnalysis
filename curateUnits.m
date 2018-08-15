function [outputArg1,outputArg2] = curateUnits(rootName, sortedData)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
load([rootName,'_all_channels.mat'],'sortedData')
load([rootName,'analysis.mat'],'Conditions','Triggers')

Spikes={};
Names={};

for i=1:size(sortedData,1)
    Spikes{i}=cell2mat(sortedData(i,2));
    Names{i}=sortedData(i,1);
end

mech=Triggers.whisker;
light=Triggers.light;

lenSpks = length(Spikes);
consIdxs = true(1,lenSpks);
% Some clusters were eliminated after an induvidual inspection of their
% PSTHs.
bads = [];
consIdxs(bads) = false;
crscor = zeros(lenSpks,lenSpks,3);    % Square matrix with Ncl x Ncl
for ccl = 1:lenSpks
    xccl = ccl+1;
    % Cross correlation avoiding the autocorrelations and the 'bads'.
    while consIdxs(ccl) && xccl <= lenSpks
        % auxXCOR = xcorr(Spikes{ccl},Spikes{xccl},'coeff');
        if consIdxs(xccl)
            clstrInfo =...
                ['Cluster ',num2str(ccl),' against cluster ',num2str(xccl)];
            disp(clstrInfo)
            % Signal reconstruction
            auxSignal1 = false(1,length(mech));
            auxSignal2 = false(1,length(mech));
            % Spikes assignment
            auxSignal1(round(fs*Spikes{ccl})) = true;
            auxSignal2(round(fs*Spikes{xccl})) = true;
            % Selecting a subsampled version of the signals given a maximum
            % lag.
            MxLag = 25; % Seconds
            rIdx = randi([round(fs*MxLag),length(mech)-round(fs*MxLag)],1);
            disp(['Random time selected: ',num2str(rIdx/fs),' seconds'])
            rWindIdx = rIdx-round(fs*MxLag):rIdx+round(fs*MxLag);
            auxSignal1 = auxSignal1(rWindIdx);
            auxSignal2 = auxSignal2(rWindIdx);
            
            % Distance matrix
            dfMtx = log(distmatrix(Spikes{ccl}',Spikes{xccl}')+1);
            lnIdx = dfMtx < log(16);
            [y,x]=find(lnIdx);
            [mdl,yhat,r2] = fit_poly(x,y,1);
            eqInfo = ['y = ',num2str(mdl(1)),'x ',num2str(mdl(2))];
            display(eqInfo)
%             figure('Name',clstrInfo);
%             imagesc(dfMtx);hold on;plot(x,yhat,'LineStyle','--',...
%                 'LineWidth',3,'Color',[1,0,0]);title([eqInfo,' ',num2str(r2)])
            crscor(ccl,xccl,1:2) = mdl;
            crscor(ccl,xccl,3) = r2;
        end
        xccl = xccl + 1;
    end
end


end

