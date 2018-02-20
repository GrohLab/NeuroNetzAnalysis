function [ dist ] = distmatrix( P1,P2,n )
%DISTMATRIX Distance matrix between each point of vectors P1 and P2.
%   n is the norm from 1 to a reasonable number (not inf included).
if nargin == 2
    n = 2;
elseif nargin > 3
    error('This function just accepts 2 input vectors, and the norm order.')
end
[points,dimension]=size(P1);
[points2,dimension2]=size(P2);
if dimension == dimension2
    dist = zeros(points,points2);
    for pf=1:points
        for pm=1:points2
            dist(pf,pm) = sum((P2(pm,:)-P1(pf,:)).^n)^(1/n);
        end
    end
else
    error(['Vector 2 has ',num2str(dimension2),...
        ' dimensions, and vector 1 has ',num2str(dimension),...
        ' dimensions'])
end
end

