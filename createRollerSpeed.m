function [iOk, vf, rollTx, fr, Texp] = createRollerSpeed(behDir)
iOk = false;
fnOpts = {"UniformOutput", false};
if ~exist(behDir, "dir")
    fprintf(1, "%s doesn't exist!\n", behDir)
    return
end
eFiles = dir(fullfile(behDir, "Roller_position*.csv"));
vFiles = dir(fullfile(behDir, "*.avi"));
Ne = numel(eFiles);
Nv = numel(vFiles);
if Ne ~= Nv
    fprintf(1, "# encoder/arduino files different from video files!\n")
    fprintf(1, "Encoder %d Video %d\n", Ne, Nv)
    return
end
%% Read and create roller 
[dt, dateFormStr] = getDates(eFiles, 'Roller_position');
crtName = @(x) fullfile(behDir, x);
getName = @(x, y) crtName(x(y).name);
dy = dt(1); dy.Format = 'yyyy-MM-dd'; dt.Format = '''T''HH_mm_ss';
auxStr = [];
for cdt = 1:numel(dt)-1
    auxStr = [auxStr, sprintf('%s+', dt(cdt))];
end
auxStr = [auxStr, sprintf('%s', dt(end))];
rsName = crtName(sprintf('RollerSpeed%s%s.mat', dy, auxStr));
if exist(rsName,"file")
    fprintf(1, "File exists! No new file created!\n")
    return
end
fr = arrayfun(@(x) VideoReader(getName(vFiles,x)), 1:Nv, fnOpts{:});
fr = cellfun(@(x) x.FrameRate, fr);
dt.Format = dateFormStr;
efName = string(arrayfun(@(x) getName(eFiles, x), 1:Ne, fnOpts{:}));
rp = arrayfun(@(x) readRollerPositionsFile(x), efName, fnOpts{:});
[vf, rollTx] = arrayfun(@(x) getRollerSpeed(rp{x}, fr(x)), 1:Ne, fnOpts{:});
vf = cellfun(@(x) x(:), vf, fnOpts{:}); vf = cat(1, vf{:});
Ns = cellfun(@(x) numel(x)-1, rollTx); 
Texp = cellfun(@(x) diff(x([1,end-1])), rollTx);
if ~std(fr)
    fr = mean(fr);
end
save(rsName, "vf","Texp","Ns","rp","rollTx","fr"); iOk = true;
end