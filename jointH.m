function hxy = jointH(im1, im2)

[row,col] = size(im1);
hxy = ones(round(max(im1(:)))+1,round(max(im2(:)))+1);
for x = 1:row
    for y = 1:col
        hxy(round(im1(x,y))+1,round(im2(x,y))+1) =...
            hxy(round(im1(x,y)+1),round(im2(x,y))+1) + 1;
    end
end
hxy = hxy/sum(hxy(:));
end

