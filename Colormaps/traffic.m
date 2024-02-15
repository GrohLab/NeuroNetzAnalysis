function [clrMap] = traffic(Nlvls)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;

checkLevels =...
    @(x) any([isnumeric(x), ~(round(x) - x), x > 0]);

addParameter(p, 'Nlvls', 256, checkLevels)

if exist("Nlvls", "var")
    parse(p, Nlvls);
else
    parse(p);
end

Nlvls = p.Results.Nlvls;

%% Auxiliary variables
getM = @(p1,p2) (p2(:,2) - p1(:,2)) ./ (p2(:,1) - p1(:,1));
getB = @(p, m) p(:,2) - m.*p(:,1);
%% Channels

hlf = Nlvls/2;
% First point
fstPts = [ones(3,1), [1;0;0]];
% Middle point
mdlPts = cat(2, hlf*ones(3,1), ones(3,1));
% Last point
lstPts = [repmat(Nlvls,3,1), [0;1;0]];

m1 = getM(fstPts, mdlPts);
m2 = getM(mdlPts, lstPts);
b1 = getB(fstPts, m1);
b2 = getB(mdlPts, m2);
clrDom = 1:Nlvls;

clrMap = cat(2, m1*clrDom(clrDom < hlf) + b1,...
    m2*clrDom(clrDom >= hlf) + b2);
clrMap = clrMap';

end

