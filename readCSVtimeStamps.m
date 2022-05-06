function csvMat = readCSVtimeStamps(fPath)
clean4Date = @(x) extractBefore(extractAfter(x,','), '+');
clean4LV = @(x) extractBefore(x,',');
dateFormStr = 'uuuu-MM-dd''T''HH:mm:ss.SSSS';
dtOpts = {'InputFormat', dateFormStr};
if exist(fPath, 'file')
    fID = fopen(fPath, 'r');lns = textscan(fID,'%s',Inf);[~] = fclose(fID);
    trigTime = datetime(clean4Date(lns{:}), dtOpts{:});
    trigTime = second(trigTime, "secondofday");
    Val = str2double(clean4LV(lns{:}));
    csvMat = [Val, trigTime];
end
end