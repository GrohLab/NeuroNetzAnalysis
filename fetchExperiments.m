function [fFlag, expFileNames] =...
    fetchExperiments(cellTypes, RecDB, anaMatDir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Nx = size(RecDB,1);
fFlag = false(Nx,1);
for cct = 1:numel(cellTypes)
    fFlag(RecDB.PhysioNucleus == cellTypes(cct)) = true;
end
expNames = RecDB.Properties.RowNames(fFlag);
expFileNames = cellfun(@fullfile,...
            repmat({anaMatDir},sum(fFlag),1),...
            cellfun(@strcat,expNames,...
            repmat({'analysis.mat'},sum(fFlag),1), 'UniformOutput',false),...
            'UniformOutput',false);
end

