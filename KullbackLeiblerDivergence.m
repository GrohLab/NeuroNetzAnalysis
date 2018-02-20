function kbd = KullbackLeiblerDivergence(P1,P2)
% KULLBACKLEIBLERDIVERGENCE returns the similarity measure between two
% probability distributions
if P1.Integrate(P1.Domain.Min,P1.Domain.Max) ~= 1.0
    P1.PDF = P1.Normalize();
end
if length(P1.PDF) == length(P2.PDF)
    N = length(P1.PDF);
    if sum(P2.PDF) ~= 1.0
        P2.PDF = P2.Normalize();
    end
elseif length(P1.PDF) > length(P2.PDF)
    warning(['The probability distribution 1 has a greater domain.\n',...
        'Attempting resampling...\n'])
    P1.PDF = P1.cutAndDownsample(P2.Domain.Min,P2.Domain.Max,length(P2.PDF));
    P1.PDF = P1.Normalize();
else
    warning(['The probability distribution 2 has a greater domain.\n',...
        'Attempting resampling...\n'])
    P2.PDF = P2.cutAndDownsample(P1.Domain.Min,P1.Domain.Max,length(P1.PDF));
    P2.PDF = P2.Normalize();
end
kbd = dot(P1.PDF,log(P1.PDF ./ P2.PDF))/log(N);
end