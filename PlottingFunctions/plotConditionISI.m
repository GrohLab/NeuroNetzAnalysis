function [genISI,h_norm,lxax] = plotConditionISI(spT, allEvents, IDe, fs, expName)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
X_AXIS_RESOLUTION = 128;
[Nn, Nsp] = size(spT);

if iscell(allEvents)
    Ne = numel(allEvents);
end
genISI = diff(spT,1,2);
genISI = (1e3*genISI)/fs;
lxax = exp(linspace(-0.7,9,X_AXIS_RESOLUTION));
invFig = figure('Visible','off');
h = histogram(log(genISI),log([0,lxax]),...
    'Visible','off','Parent',invFig);
h_norm = h.Values/sum(h.Values);

close(invFig)
figure('Name',expName,'Color',[1,1,1]);
ax(1) = subplot(1,2,1);semilogx(lxax,h_norm,...
    'LineWidth',3,'Color',[0,0,0],'DisplayName','p(\Deltat))');
title('Histogram ISI')
xlabel('Time [ms]');ylabel('p(\Deltat)');grid('minor')
ax(2) = subplot(1,2,2);semilogx(lxax,cumsum(h_norm),...
    'LineWidth',3,'Color',[0,0,0],'DisplayName','p(\Deltat))');
title('Cumulative histogram ISI')
xlabel('Time [ms]');ylabel('p(\Deltat)');grid('minor');
evntIdx = find(cellfun(@sum,allEvents) > 0);
emptyEvents = sum(cellfun(@sum,allEvents) == 0);
trueNe = Ne - emptyEvents;
h_cond = zeros(trueNe,X_AXIS_RESOLUTION);
clrMap = bone(trueNe);
for ce = 1:trueNe
    auxSpkLabel = allEvents{evntIdx(ce)}(spT);
    if ~strcmp(IDe{evntIdx(ce)},'exclude') && ~strcmp(IDe{evntIdx(ce)},'grooming')
        auxSpkLabel = allEvents{evntIdx(ce)}(spT);
        auxFig = figure('Visible','off');
        h = histogram(log(genISI(auxSpkLabel(1:end-1))),log([0,lxax]),...
            'Visible','off','Parent',auxFig);
        h_cond(ce,:) = h.Values/sum(h.Values);
        
        close(auxFig)
        if ce == 1
            hold(ax(1),'on')
            hold(ax(2),'on')
        end
        semilogx(ax(1),lxax,h_cond(ce,:),'DisplayName',IDe{evntIdx(ce)},...
            'Color',clrMap(ce,:))
        semilogx(ax(2),lxax,cumsum(h_cond(ce,:)),'DisplayName',IDe{evntIdx(ce)},...
            'Color',clrMap(ce,:))
        
    else
        fprintf('Excluding %d spikes due to a less controled environment.\n',...
            sum(auxSpkLabel))
    end
end

legend(ax(1),'show')
legend(ax(2),'show')

end

% h_norm = zeros(Nn,X_AXIS_RESOLUTION,'single');
% for cn = 1:Nn
%     h = histogram(log(genISI),X_AXIS_RESOLUTION,...
%         'Visible','off','Parent',invFig);
%     h_norm(cn,:) = h.Values/sum(h.Values);
%     subplot(1,2,1);semilogx(lxax,h_norm);title('Histogram ISI')
%     xlabel('Time [s]');ylabel('p(\Deltat)');grid('minor')
%     subplot(1,2,2);semilogx(lxax,cumsum(h_norm));title('Cumulative histogram ISI')
%     xlabel('Time [s]');ylabel('p(\Deltat)');grid('minor')
% end
