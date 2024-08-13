function dlcTable = readDLCData(dlcPath)
%READDLCDATA reads the .csv file from DeepLabCuts and creates a table with
%all body parts and its estimated x and y coordinates, and its likelihood
%at that time point.
fnOpts = {'UniformOutput', false};
fprintf(1, 'Reading file... ');
dlcData = importdata(dlcPath, ',', 3);
fprintf(1, 'Ready!\n')
bdyPrts = unique( strsplit( dlcData.textdata{2,1}, ',' ), 'stable' );
bdyPrts(1) = [];
% Packing x, y, and likelihood in a cell
xylData = arrayfun(@(x) dlcData.data(:,x:x+2),...
    (2:3:size(dlcData.data,2)-2), fnOpts{:});
dlcTable = table(xylData{:}, 'VariableNames', bdyPrts);
end