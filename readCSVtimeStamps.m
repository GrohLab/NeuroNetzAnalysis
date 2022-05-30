function csvMat = readCSVtimeStamps(fPath)
clean4Date = @(x) extractBefore(extractAfter(x,','), '+');
clean4LV = @(x) extractBefore(x,',');
dateFormStr = 'uuuu-MM-dd''T''HH:mm:ss.SSSS';
dtOpts = {'InputFormat', dateFormStr};
if exist(fPath, 'file')
    fID = fopen(fPath, 'r');lns = textscan(fID,'%s',Inf);[~] = fclose(fID);
    trigTime = datetime(clean4Date(lns{:}), dtOpts{:});
    trigTime = second(trigTime, "secondofday");
    Val = arrayfun(@(x) str2num(x,"Evaluation","restricted"), ...
        lower(clean4LV(string(lns{:})))); %#ok<ST2NM> 
    csvMat = [Val, trigTime];
end
end