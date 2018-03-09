function [graphObjs] = plotComparisonBars(cond1,cond2,str1,str2)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

graphObjs(1) = plot(ones(1,numel(cond1)),cond1,'LineStyle','none','Marker','.');
hold('on')
graphObjs(2) = plot(2*ones(1,numel(cond2)),cond2,'LineStyle','none','Marker','.');
errorbar(0.9,mean(cond1),std(cond1),'Color',graphObjs(1).Color,'Marker','o')
errorbar(2.1,mean(cond2),std(cond2),'Color',graphObjs(2).Color,'Marker','o')
xticks([1,2]);xticklabels({str1,str2});
end

