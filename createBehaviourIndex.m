function [pAreas, vertices, vFig] = createBehaviourIndex(behRes,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
p = inputParser;
BSfields = {'ConditionName', 'Results'};
fnOpts = {'UniformOutput', false};

validateBS = @(x) isstruct(x) & all(isfield(x, BSfields));
addRequired(p, 'behRes', validateBS)

parse(p, behRes, varargin{:})

behRes = p.Results.behRes;

%%

consCondNames = string({behRes.ConditionName});
Nccond = numel( consCondNames ); clrMp = lines(Nccond);
bsNames = {behRes(1).Results.BehSigName}; 
mvPr = arrayfun(@(br) [br.Results.MovProbability], behRes, fnOpts{:});
mvAm = arrayfun(@(br) [br.Results.AmplitudeIndex], behRes, fnOpts{:});
mv = {mvPr, mvAm}; figNames = ["Trial proportion", "Amplitude index"];
pAreas = zeros(Nccond, 2); vFig = gobjects(2, 1);
for cmv = 1:numel(mv)
    [bAxis, vertices] = prepareRadialData( behRes, mv{cmv} );
    
    [vFig(cmv), pd] = prepareRedialGraphics( bAxis, bsNames, figNames(cmv) );
    % plot data (pd)
    pAreas(:,cmv) = pd(vertices, clrMp, consCondNames);
end

end

function [bAxis, vertices, Nccond] = prepareRadialData(behRes, medi)
fnOpts = {'UniformOutput', false};
Nccond = numel(behRes);
mvPr = cat(1, medi{:});
Nb = size(mvPr, 2); radAxis = (0:Nb-1)*(2*pi/Nb); 
bAxis = exp(1i*radAxis(:));
vertices = arrayfun(@(c) mvPr(c,:).*(bAxis.'), 1:Nccond, fnOpts{:});
vertices = cat(1, vertices{:});
end

function [vFig, fpd] = prepareRedialGraphics(bAxis, bsNames, figName)
axOpts = {'Box', 'off', 'Color', 'none'};
lgOpts = [axOpts(:)', {'Location'}, {'best'}];
txOpts = {'HorizontalAlignment','center', 'VerticalAlignment', 'middle', ...
    'Rotation'};
vFig = figure('Name', figName, 'Color', 'w'); 
ax = axes('Parent', vFig, 'NextPlot', 'add', axOpts{:});
grayLvl = 0.75; gridThick = 0.1;
% R-grid
line([zeros(size(bAxis)), real(bAxis)]', ...
    [zeros(size(bAxis)), imag(bAxis)]', 'LineWidth', gridThick, ...
    'Color', grayLvl*ones(1,3));
% Theta-grid
tcks = (1:4)/4; th = 0:pi/50:2*pi; 
x_circ = tcks(:) .* cos(th); y_circ = tcks(:) .* sin(th);
line( x_circ', y_circ', 'LineWidth', gridThick, 'Color', grayLvl*ones(1,3) )
arrayfun(@(c,t,a) text(real(c)*1.1, imag(c)*1.1, t, txOpts{:}, a ), ...
    bAxis, string( bsNames(:) ), (180*angle( bAxis )/pi) - 90 );

% arrayfun(@(r) text(real(bAxis(r))*tcks, imag(bAxis(r))*tcks,sprintf("%c",743), ...
%     txOpts{:}, rad2deg(angle(bAxis(r)))), 1:size(bAxis,1));
text(tcks, zeros(numel(tcks),1), string(tcks), txOpts{3}, 'top' )
set(ax, axOpts{:})
fpd = @plotData;
    function pAreas = plotData(vertices, clrMp, consCondNames)
        [Nccond, Nb] = size(vertices);
        pObj = arrayfun(@(r) patch( real( vertices(r,:) ), ...
            imag( vertices(r,:) ), ones( 1, Nb ), ...
            'FaceAlpha', 0.1, 'EdgeColor', clrMp(r,:), ...
            'FaceColor', clrMp(r,:), 'LineStyle', '-.'), 1:Nccond);
        pShape = arrayfun(@(p) polyshape(p.Vertices(:,1), p.Vertices(:,2)), pObj);
        pAreas = area(pShape); set(ax, axOpts{:}, 'Visible', 'off')
        axis(ax, [-1,1,-1,1], 'square'); 
        lgObj = legend(pObj, consCondNames + " " + string(round(pAreas,3))); 
        set(lgObj, lgOpts{:})
    end
end


