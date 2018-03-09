function  [vctrs,measures] = eigenAnalysis(analiticalSignals,idxs,normFlag)
%EIGENANALYSIS Summary of this function goes here
%   Detailed explanation goes here
[Ns,samples] = size(analiticalSignals);
if Ns > samples
    analiticalSignals = analiticalSignals';
    analiticalSignals = conj(analiticalSignals);    
end
if nargin < 3
    normFlag = false;
    if nargin < 2
        idxs = true(1,samples);
    end
else
    if isempty(idxs)
        idxs = true(1,samples);
    end
end
sigCS = zeros(2,2,Ns);
vctrs = zeros(2,2,Ns);
measures = zeros(3,Ns);
for cs = 1:Ns
    % Normalizing each signal to the infinity norm.
    if normFlag
        currentSignal = analiticalSignals(cs,:)/...
            norm(analiticalSignals,'inf');
    else
        currentSignal = analiticalSignals(cs,:);
    end
    sigCS(:,:,cs) = cov(real(currentSignal(idxs)),...
        imag(currentSignal(idxs)));
    [vctrs(:,:,cs),L] = eig(sigCS(:,:,cs));
    vctrs(:,:,cs) = vctrs(:,:,cs)*L;
    L = diag(L);
    % Measuring the spread of the points
    measures(:,cs) = [...
        max(L)\min(L);... 
        hypot(L(1),L(2));...
        atan2(vctrs(2,2,cs),vctrs(1,2,cs));...
        ];
end

