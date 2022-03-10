function iOk = createRollerSpeed(behDir)
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
crtName = @(x) fullfile(behDir, x);
getName = @(x, y) crtName(x(y).name);
fr = arrayfun(@(x) VideoReader(getName(vFiles,x)), 1:Nv, fnOpts{:});
fr = cellfun(@(x) x.FrameRate, fr);
[dt, dateFormStr] = getDates(eFiles, 'Roller_position');
dt.Format = dateFormStr;
efName = string(arrayfun(@(x) getName(eFiles, x), 1:Ne, fnOpts{:}));
rp = arrayfun(@(x) readRollerPositionsFile(x), efName, fnOpts{:});
[vf, rollTx] = arrayfun(@(x) getRollerSpeed(rp{x}, fr(x)), 1:Ne, fnOpts{:});
end