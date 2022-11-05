function linTrans = scaleOnLine(data,y1,y2)
x1 = min(data(:));
x2 = max(data(:));
m = (y2 - y1)/(x2 - x1);
b = y2 - m*x2;
linTrans = m.*data + b;
end