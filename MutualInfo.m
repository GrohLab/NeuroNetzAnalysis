function [ MI ] = MutualInfo( im1, im2 )

hxy = jointH(im1, im2);
[x,y]=size(hxy);
MI = (marginalEnt(hxy,2) + marginalEnt(hxy,1) - ...
    sum(sum(hxy.*log10(hxy))))/log2(x*y);

end