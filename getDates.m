function [fDates, dateFormStr] = getDates(fileNames, baseName, endng)
% Function to extract the date from the file name.
% Extract machine readable 19-December-2021 15:32:12 from
% Roller_position2021-12-19T15_32_12.csv, for example. "Roller_position"
% would be the baseName, and an array with structures from dir function.
dateFormStr = 'yyyy-MM-dd''T''HH_mm_ss';
fnOpts = {'UniformOutput', false};
if ~exist("endng", "var")
    endng = ".";
end
getDate = @(x, y) extractBefore(extractAfter(x, y), endng);
switch class(fileNames)
    case 'struct'
        fDates = arrayfun(@(x) getDate(x.name, baseName), fileNames,...
            fnOpts{:});
    case {'char', 'string'}
        fDates = arrayfun(@(x) getDate(x, baseName), fileNames,...
            fnOpts{:});
    otherwise
        error('Unrecognised input format! Must be a string or dir() output!')
end
try
    fDates = cellfun(@(x) datetime(x, 'InputFormat', dateFormStr),...
        fDates);
    fDates.Format = dateFormStr;
catch ME
    disp(ME)
    fprintf(1, "Returning string, not date!\n");
    fDates = [fDates{:}];
end
end