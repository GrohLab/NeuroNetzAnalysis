function [trigTime, trigName] = getCSVTriggers(fPath)
clean4Date = @(x) extractBefore(extractAfter(x,','), '+');
clean4LV = @(x) extractBefore(x,',');
dateFormStr = 'uuuu-MM-dd''T''HH:mm:ss.SSSS';
dtOpts = {'InputFormat', dateFormStr};
fnOpts = {'UniformOutput', false};
if exist(fPath, 'file')
    fID = fopen(fPath, 'r');lns = textscan(fID,'%s',Inf);[~] = fclose(fID);
    trigTime = datetime(clean4Date(lns{:}), dtOpts{:});
    trigTime = second(trigTime, "secondofday");
    udVal = arrayfun(@str2num, lower(string(clean4LV(lns{:}))));
    if sum(udVal) == sum(~udVal)
        trigTime = cat(2, trigTime(udVal), trigTime(~udVal));
    else
        trigTime = arrayfun(@(x) trigTime(xor(udVal, x)), [true false], fnOpts{:});
    end
    if ~udVal(1)
        trigTime = flip(trigTime,2);
    end
    [~, trigName, ~] = fileparts(fPath); trigName = extract(trigName,1);
end
end