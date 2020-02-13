function WhiteData = dataWhitening(data)
dataZ = data - mean(data,1);
sig = cov(data);
[U,S] = svd(sig);
Sinv = diag(1./sqrt(diag(S)));
WhiteData = (Sinv * U' * dataZ')';

end