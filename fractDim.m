function [FD, r2 ] = fractDim( signal, K )
%FRACTDIM This function returns the fractal dimension of the input SIGNAL
%using the Higuchi's algorithm.
%   FD = fractDim(SIGNAL, K) returns the fractal dimension FD considering
%   the starting point m until reaching the $$((N-m)/K)*K$$ sample. Where N
%   is the total number of samples in the signal.
mL = zeros(K,1);
N=length(signal);
for k=1:K
    mL(k) = meanLength(signal,k,N);
end

lL = log(mL);
lK = log(1:K);

[linFit,FD]=detrend_profile(1,lK,lL');
r2 = fitGoodness(linFit,lL');

end

function mL = meanLength(signal,K,N)
m = 1;
idxs = m:K:(m + round((N-m)/K)*K);
x_mk = zeros(K,length(idxs));
n_orm = zeros(K,1);
L_m = zeros(K,1);
for m=1:K
    if (m + round((N-m)/K)*K)>N
        idxs = m:K:N;
    else
        idxs = m:K:(m + round((N-m)/K)*K);
    end
    x_mk(m,1:length(idxs)) = signal(idxs);
    n_orm(m) = ((N-1)/((round(N-m)/K)*K))/K;
    L_m(m) = n_orm(m) * sum(abs(diff(x_mk(m,:))));
end
mL = mean(L_m);
end