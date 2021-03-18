function [iOk] = exportSpks2Txt(spkSubs, saveDir, fs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
milliFact = 1e3/fs;
spkTms = cellfun(@(x) x*milliFact, spkSubs, 'UniformOutput', false);
formatSpec = '%f\n'; iOk = false;
if ~exist(saveDir,'dir')
    fprintf(1, 'Given directory doesn''t exist!\n');
    fprintf(1, 'Trying to create it... ');
    [mkOk, ~, ~] = mkdir(saveDir);
    if mkOk
        fprintf(1, ' success!\n');
    else
        fprintf(1, ' failed. Returning...\n');
        return
    end
end
for a = 1:length(spkTms)
    ftxtName = fullfile(saveDir, "cell"+(a-1)+".txt");
    fileID = fopen(ftxtName, 'w');
    if fileID >= 3
        fprintf(fileID, formatSpec, spkTms{a}); fclose(fileID);
    else
        fprintf(1, 'Not able to create/open %s file!\n', ftxtName)
    end
end
iOk = true;
end

