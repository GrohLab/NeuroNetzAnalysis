function pDir = getParentDir(dir,p)
auxDirs = strsplit(dir,'\');
pieces = length(auxDirs);
if isempty(auxDirs{pieces})
    auxDirs(pieces) = [];
    pieces = pieces - 1;
end
pDir = [];
for d = 1:pieces-p
    pDir = [pDir,auxDirs{d},'\'];
end
end
