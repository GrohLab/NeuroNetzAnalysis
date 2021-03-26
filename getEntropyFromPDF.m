function e = getEntropyFromPDF(P)
N = length(P);
if sum(P)~=1
    P = P/sum(P);
end
P(P==0) = range(P)*eps;

e = -dot(P,log(P))/log(N);
end