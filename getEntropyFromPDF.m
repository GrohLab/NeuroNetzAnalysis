function e = getEntropyFromPDF(P)
N = length(P);
if sum(P)~=1
    P = P/sum(P);
end
e = -dot(P,log(P))/log(N);
end