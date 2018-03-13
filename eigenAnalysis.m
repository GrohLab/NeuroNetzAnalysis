function  [vctrs,measures] = eigenAnalysis(analiticalSignals,idxs,normFlag)
%EIGENANALYSIS Computes the eigenvectors and values for the given
%analitical signal.
%   [vectrs, measures] = eigenAnalysis(analiticalSignal, indexes, normFlag)
%   The function accepts 1, 2, or 3 input arguments. The analitical signal
%   is the signal to be considered and processed, the indexes are to select
%   specific samples in the signal e. g. whisking spikes, triggered evoked
%   spikes, etc., and the normFlag is recommended to be set to True in
%   phase-locking value computations. Otherwise, the normalization step can
%   be skipped.

[Ns,samples] = size(analiticalSignals);
if Ns > samples
    analiticalSignals = analiticalSignals';
    analiticalSignals = conj(analiticalSignals);    
end
if nargin < 3
    normFlag = true;
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
        currentSignal = analiticalSignals(cs,:)./...
            abs(analiticalSignals(cs,:));
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

