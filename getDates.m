function [fDates, dateFormStr] = getDates(fileNames, baseName)
% Function to extract the date from the file name.
% Extract machine readable 19-December-2021 15:32:12 from
% Roller_position2021-12-19T15_32_12.csv, for example. "Roller_position"
% would be the baseName, and an array with structures from dir function.
dateFormStr = 'yyyy-MM-dd''T''HH_mm_ss';
fnOpts = {"UniformOutput", false};
getDate = @(x, y) extractBefore(extractAfter(x, y), ".");
fDates = arrayfun(@(x) getDate(x.name, baseName), fileNames,...
    fnOpts{:});
fDates = cellfun(@(x) datetime(x, 'InputFormat', dateFormStr),...
    fDates);
end