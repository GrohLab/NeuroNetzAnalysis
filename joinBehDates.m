function [outName] = joinBehDates(inFiles, inFilePttrn, ...
    outFilePttrn, endng)
%JOINBEHDATES creates the name for the joint file using the date and time
%of the input files.
%   Detailed explanation goes here
Nif = numel(inFiles);
fixpath = @(x) insertAfter(x,"\","\");
if ~exist("endng", "var")
    endng = ".";
end
[dt, dateFormStr] = getDates(inFiles, inFilePttrn, endng);
if isdatetime(dt)
    dy = dt(1); dy.Format = 'yyyy-MM-dd'; 
    dt.Format = '''T''HH_mm_ss';
    auxStr = "";
    for cdt = 1:Nif-1
        auxStr = auxStr + sprintf("%s+", dt(cdt));
    end
    auxStr = auxStr + sprintf("%s", dt(end));
    outName = sprintf(fixpath(outFilePttrn), dy, auxStr);
    if exist(outName,"file")
        fprintf(1, "File exists!\n")
    end
else
    fprintf(1, "Debugging needed!\n")
    dbstop in joinBehDates at 23
end
