function [pAreas, vertices, vFig] = createBehaviourIndex(behRes,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
p = inputParser;
BSfields = {'ConditionName', 'Results'};

validateBS = @(x) isstruct(x) & all(isfield(x, BSfields));
addRequired(p, 'behRes', validateBS)

parse(p, behRes, varargin{:})

behRes = p.Results.behRes;

%%

consCondNames = string({behRes.ConditionName});
bsNames = {behRes(1).Results.BehSigName}; 
[bAxis, vertices, Nccond, axTtlPos] = prepareRadialData(behRes);
clrMp = lines(Nccond);

[vFig, pd] = prepareRedialGraphics(bAxis, bsNames, axTtlPos);
% plot data (pd)
pAreas = pd(vertices, clrMp, consCondNames);

end

function [bAxis, vertices, Nccond, btmx] = prepareRadialData(behRes)
fnOpts = {'UniformOutput', false};
Nccond = numel(behRes);
mvPr = arrayfun(@(br) [br.Results.MovProbability], behRes, fnOpts{:});
mvPr = cat(1, mvPr{:});
Nb = size(mvPr, 2); radAxis = (0:Nb-1)*(2*pi/Nb); 
bAxis = exp(1i*radAxis(:)); btmx = 1.15*max(mvPr,[],'all');
vertices = arrayfun(@(c) mvPr(c,:).*transp(bAxis), 1:Nccond, fnOpts{:});
vertices = cat(1, vertices{:});
end

function [vFig, fpd] = prepareRedialGraphics(bAxis, bsNames, btmx)
axOpts = {'Box', 'off', 'Color', 'none'};
lgOpts = [axOpts(:)', {'Location'}, {'best'}];
vFig = figure('Name', 'Behaviour dimensions', 'Color', 'w'); 
ax = axes('Parent', vFig, 'NextPlot', 'add', axOpts{:});
plot([zeros(size(bAxis)), real(bAxis)]', ...
    [zeros(size(bAxis)), imag(bAxis)]', 'LineWidth', 0.3, ...
    'Color', 0.3*ones(1,3), 'LineStyle', ':'); 
text(real(bAxis)*btmx*1.4, imag(bAxis)*btmx*1.4, string(bsNames), ...
    'HorizontalAlignment','center')
tcks = round(btmx,1)*(1:4)/4;
arrayfun(@(r) text(real(bAxis(r))*tcks, imag(bAxis(r))*tcks,'|', ...
    'Rotation', rad2deg(angle(bAxis(r))), ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle'), 1:size(bAxis,1)); 
text(tcks, zeros(numel(bAxis),1), string(tcks), 'HorizontalAlignment','center', ...
    'VerticalAlignment','top')
set(ax, axOpts{:})
fpd = @plotData;
    function pAreas = plotData(vertices, clrMp, consCondNames)
        [Nccond, Nb] = size(vertices);
        pObj = arrayfun(@(r) patch(real(vertices(r,:)), imag(vertices(r,:)), ones(1,Nb), ...
            'facealpha', 0.3, 'EdgeColor', 'none', 'FaceColor', clrMp(r,:)), ...
            1:Nccond);
        pShape = arrayfun(@(p) polyshape(p.Vertices(:,1), p.Vertices(:,2)), pObj);
        pAreas = area(pShape); set(ax, axOpts{:}, 'Visible', 'off')
        axis(ax, (2*btmx*(btmx<0.5) + btmx>0.5)*[-1,1,-1,1], 'square'); 
        lgObj = legend(pObj, consCondNames + " " + string(round(pAreas,3))); 
        set(lgObj, lgOpts{:})
    end
end


