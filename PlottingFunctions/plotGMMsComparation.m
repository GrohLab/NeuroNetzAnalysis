function [fig,p_xhat,KLd] = plotGMMsComparation(parameters,Xdomain,IDe)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
normPDF = @(x) x./sum(x);
isaline = @(x) isa(x,'matlab.graphics.chart.primitive.Line');

M = size(parameters,1);
C = size(parameters,3);
pik = zeros(M,length(Xdomain),C);
for cc = 1:C
    if all(all(isnan(parameters(:,:,cc))))
        continue
    end
    for k=1:M
        if all(isnan(parameters(k,:,cc)))
            continue
        end
        pik(k,:,cc)=evalgauss(Xdomain,parameters(k,2,cc),parameters(k,3,cc));
        pik(k,:,cc)=pik(k,:,cc).*parameters(k,1,cc);
    end
end
p_xhat = squeeze(sum(pik,1))';
conComb = combnk(1:C,2);
Ncmp = size(conComb,1);
KLd = zeros(Ncmp,3);

for cpdf = 1:Ncmp
    KLd(cpdf,1:2) = conComb(cpdf,:);
    if ~sum(p_xhat(conComb(cpdf,1),:)) || ~sum(p_xhat(conComb(cpdf,2),:))
        continue
    end
    KLd(cpdf,3) = KullbackLeiblerDivergence(...
        p_xhat(conComb(cpdf,1),:),...
        p_xhat(conComb(cpdf,2),:));
end
cmap = lines(C);
fig = figure('Name','GMMs for the given conditions','Color',[1,1,1]);
lns = gobjects(C,1);
for cpdf = 1:C
    if sum(p_xhat(cpdf,:)) == 0
        continue
    end
    p_xhat(cpdf,:) = normPDF(p_xhat(cpdf,:));
    if cpdf > 1
        ax = fig.Children;
        set(ax,'NextPlot','add')
    end
    yyaxis('left')
    if exist('IDe','var') && ~isempty(IDe)
        lns(cpdf) = plot(Xdomain,p_xhat(cpdf,:),'LineStyle','-',...
            'Color',cmap(cpdf,:),...
            'DisplayName',IDe{cpdf},...
            'Marker','none');
    else
        lns(cpdf) = plot(Xdomain,p_xhat(cpdf,:),'LineStyle','-',...
            'Color',cmap(cpdf,:),...
            'DisplayName',sprintf('Condition %d',...
            cpdf),...
            'Marker','none');
    end
    yyaxis('right')
    plot(Xdomain,cumsum(p_xhat(cpdf,:)),'Color',cmap(cpdf,:),...
        'LineStyle','--', 'Marker','none')
end
lFlag = arrayfun(isaline,lns);
title('PDFs comparison');box off
legend(lns(lFlag))
xlabel('Time [s]');
ax.YAxis(1).Label = 'Probability';
ax.YAxis(1).Color = [0,0,0];
ax.YAxis(2).Color = [0,0,0];
ax.YAxis(2).Limits = [0,1];
ax.YAxis(2).TickValues = [0, 0.5, 1];
ax.YAxis(2).Label = 'Cumulative probability';

yyaxis('left');ylabel('Probability')
yyaxis('right');ylabel('Cumulative probability');
end

