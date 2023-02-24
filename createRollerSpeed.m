function [iOk, vf, rollTx, fr, Texp] = createRollerSpeed(behDir)
%% Auxiliary variables and functions
iOk = false;
fnOpts = {'UniformOutput', false};
if ~exist(behDir, "dir")
    fprintf(1, "%s doesn't exist!\n", behDir)
    return
end
eFiles = dir(fullfile(behDir, "Roller_position*.csv"));
vFiles = dir(fullfile(behDir, "*.avi"));
aFiles = dir(fullfile(behDir, "ArduinoTriggers*.mat"));
Ne = numel(eFiles);
Nv = numel(vFiles);
Na = numel(aFiles);
if Ne ~= Nv
    fprintf(1, "# encoder/arduino files different from video files!\n")
    fprintf(1, "Encoder %d Video %d\n", Ne, Nv)
    return
end
if Ne ~= Na
    fprintf(1, "# encoder files different from Trigger files!\n")
    fprintf(1, "Encoder %d Triggers %d\n", Ne, Na)
    fprintf(1, "Run 'readAndCorrectArdTrigs' first to create Trigger files\n")
    return
end
    function minOfSt= getArd_IntOffset()
        % Trigger correspondance between cells
        [iS, aS] = arrayfun(@(x) find(x.atNames == x.itNames), tStruct, ...
            fnOpts{:});
        % Minimum offset time between trigger signals
        ofSts = arrayfun(@(x, ao, io) arrayfun(@(at, it) ...
            abs(at{1}(1) - it{1}(1,1)),x.atTimes(ao{:}), x.itTimes(io{:})'), ...
            tStruct, aS, iS, fnOpts{:});
        minOfSt = cellfun(@min, ofSts);
    end
%% Read and create roller 
[dt, dateFormStr] = getDates(eFiles, 'Roller_position');
crtName = @(x) fullfile(behDir, x);
getName = @(x, y) crtName(x(y).name);
flfa = @(x) fullfile(x.folder, x.name);
vars2load = {'Nt', 'atTimes', 'atNames', 'itTimes', 'itNames', 'minOfSt'};
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
tStruct = arrayfun(@(x) load(flfa(x),vars2load{:}), aFiles);
fr = arrayfun(@(x) VideoReader(getName(vFiles,x)), 1:Nv, fnOpts{:});
fr = cellfun(@(x) x.FrameRate, fr);
dt.Format = dateFormStr;
efName = string(arrayfun(@(x) getName(eFiles, x), 1:Ne, fnOpts{:}));
rp = arrayfun(@(x) readRollerPositionsFile(x), efName, fnOpts{:});
[vf, rollTx] = arrayfun(@(x) getRollerSpeed(rp{x}, fr(x)), 1:Ne, ...
    fnOpts{:});
if ~isfield(tStruct, 'minOfSt')
    minOfSt = getArd_IntOffset();
else
    minOfSt = [tStruct.minOfSt];
end
rpFrmt = string(cellfun(@class, rp, fnOpts{:}));
switch rpFrmt
    case "double"
        if minOfSt > 0
            vf = arrayfun(@(x) padarray(vf{x}, [0, ...
                round(minOfSt(x)*fr(x))], 0, 'pre'), 1:Ne, fnOpts{:});
        else
            vf = arrayfun(@(x) vf{x}(round(abs(minOfSt(x))*fr(x)):end), 1:Ne, ...
                fnOpts{:});
        end
    case "table"
        ibFlags = arrayfun(@(x, y, z) (x{:} - (10 + second(y, "secondofday") ...
            - z)) < 0, rollTx, dt, minOfSt, fnOpts{:});
        vf = cellfun(@(x, y) x(~y(1:end-1)), vf, ibFlags, fnOpts{:});
    otherwise
end
Texp = arrayfun(@(x) size(vf{x}, 2)/fr(x), 1:Ne);
if isfield(tStruct, 'Nt')
    vf = arrayfun(@(x) padarray(vf{x}, [0, round(abs(Texp(x) - ...
        tStruct(x).Nt)*fr(x))], 0, 'post'), 1:Ne, fnOpts{:});
    Texp = arrayfun(@(x) length(vf{x})/fr(x), 1:Ne);
    rollTx = arrayfun(@(x) ...
        rollTx{x}(1):1/fr(x):Texp(x)+rollTx{x}(1)-1/fr(x), 1:Ne, fnOpts{:});
end
vf = cellfun(@(x) x(:), vf, fnOpts{:}); Ns = cellfun(@(x) numel(x), vf);
vf = cat(1, vf{:});
if ~std(fr)
    fr = mean(fr);
else
    fprintf(1, "Frame rates heterogeneous!\n")
end
save(rsName, "vf","Texp","Ns","rp","rollTx","fr"); iOk = true;
end