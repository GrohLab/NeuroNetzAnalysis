function [outputArg1,outputArg2] = createBehaviourIndex(behRes,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
p = inputParser;
BSfields = {'ConditionName', 'Results'};

validateBS = @(x) isstruct(x) & isfield(x, BSfields{:});
addRequired(p, 'behRes', validateBS)

parse(p, behRes, varargin{:})

behRes = p.Results.behRes;

%%
fnOpts = {'UniformOutput', false};
axOpts = {'Box', 'off', 'Color', 'none'};
lgOpts = [axOpts(:)', {'Location'}, {'best'}];

Nccond = numel(behRes);

bsNames = {behRes(1).Results.BehSigName};
mvPr = arrayfun(@(br) [br.Results.MovProbability], behRes, fnOpts{:});
mvPr = cat(1, mvPr{:}); Nb = size(mvPr, 2); radAxis = (0:Nb-1)*(2*pi/Nb); 
bAxis = exp(1i*radAxis(:)); axTtlPos = 1.15*max(mvPr,[],'all');
vertices = arrayfun(@(c) mvPr(c,:).*transp(bAxis), 1:Nccond, fnOpts{:});
vertices = cat(1, vertices{:}); clrMp = lines(Nccond);

figure; plot([zeros(size(bAxis)),real(bAxis)]', ...
    [zeros(size(bAxis)), imag(bAxis)]', 'LineWidth', 0.5, ...
    'Color', 0.25*ones(1,3)); hold on; 
text(real(bAxis)*axTtlPos,imag(bAxis)*axTtlPos, string(bsNames), ...
    'HorizontalAlignment','center')
pObj = arrayfun(@(r) patch(real(vertices(r,:)), imag(vertices(r,:)), ones(1,Nb), ...
    'facealpha', 0.3, 'EdgeColor', 'none', 'FaceColor', clrMp(r,:), ...
    'DisplayName', consCondNames(r)), 1:Nccond);
ax = gca; set(ax, axOpts{:}, 'Visible', 'off')
lgObj = legend(pObj); set(lgObj, lgOpts{:})
end

