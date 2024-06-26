function kbd = KullbackLeiblerDivergence(P1,P2)
% KULLBACKLEIBLERDIVERGENCE returns the similarity measure between two
% probability distributions
P1 = checkDist( P1 );
P2 = checkDist( P2 );
if length(P1) == length(P2)
    N = length(P1);
else
    kbd = NaN;
    fprintf('! The distributions have different resolutions (different length)!\n')
    fprintf('No divergence calculated\n')
    return;
end
aux = rdivide(P1,P2);
if any( isinf(aux) )
    kbd = dot( P1, log2(P1) - log2(P2) ) / log2(N);
else
    kbd = dot(P1,log2(aux))/log2(N);
end

end

function P = checkDist(P)
P(P == 0) = range(P)*1e-6;
P = P/sum(P);
end