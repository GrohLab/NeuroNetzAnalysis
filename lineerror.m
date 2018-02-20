function [err, ds] = lineerror(mdl,pts)
%LINEERROR returns the error and the distance to the line for the given
%points. 
%   [error, distance] = LINEERROR(model, 2-D_points);
cPt1 = pts(pts(:,1)==max(pts(:,1)),1);
cPt1=cPt1(1);
cPt2 = pts(pts(:,1)==min(pts(:,1)),1);
cPt2=cPt2(1);
cPt1(2)=mdl(1)*cPt1(1) + mdl(2);
cPt2(2)=mdl(1)*cPt2(1) + mdl(2);
r = cPt2 - cPt1;
n = [-r(2);r(1)];
n = n/norm(n,2);
d = cPt1*n;
ds = pts*n - d;
err = sum(abs(ds))/sum(abs(hypot(pts(:,1),pts(:,2))));
end
