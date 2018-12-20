%OUTPUTS
function [profile,...   Auxiliary output
    Hq,...              Hurst exponents
    tq,...              Scaling exponents
    hq,...              Singularity strength
    Dq,...              hq subset dimensionality
    Fq]...              Average over all elements of the q-order fluctuation function
    = multiFract_DFA(...%INPUTS
    signal,...          Time-series to be analyzed
    ts,...              Time axis
    scale,...           Scales to separate the signal
    q,...               q-order statistics
    m)...               m-order polynomial to detrend
    %MULTIFRACT_DFA Computation of the multifractal spectrum using $q$-orders
%   Detailed explanation goes here

% Determining the profile as in Kantelhardt et al. 2002
profile = cumsum(signal - mean(signal));
% profile = cumsum(profile - mean(profile));
Fq = zeros(length(q),length(scale));
qbs = cell(size(Fq));
bs = cell(1,length(scale));
for s = 1:length(scale)
    nseg = floor(length(profile)/scale(s));
    for v=1:nseg
        idxs = uint32((((v-1) * round(scale(s))) + 1):v*scale(s));
        [polft,~]=detrend_profile(m,ts(idxs),profile(idxs));
        bs{s}(v) = sqrt(mean((profile(idxs)-polft).^2));
    end
    % Current q (cq)
    for cq = 1:length(q)
        qbs{cq,s} = bs{s}.^q(cq);
        Fq(cq,s)=mean(qbs{cq,s}).^(1/q(cq));
    end
    % Correcting the variations for q-order zero.
    Fq(q==0,s)=exp(0.5*mean(log(bs{s}.^2)));
end

for cq = 1:length(q)
    [~,ft] = detrend_profile(1,log2(scale),log2(Fq(cq,:)));
    Hq(cq) = ft(1);
    qrL = log2(scale).*Hq(cq) + ft(2);
end
if length(q)~=1
    tq = (Hq.*q) - 1;                   % Scaling exponent
    hq = diff(tq)./(q(2)-q(1));         % Singularity strength / Hölder exp.
    Dq =(q(1:end-1).*hq)-tq(1:end-1);   % Dimension of subset characterized by hq
else
    tq = 1;
    hq = tq;
    Dq = 1;
end
end


% Kantelhardt, J. W., Zschiegner, S. A., Koscielny-Bunde, E., Havlin, S.,
% Bunde, A., & Stanley, H. E. (2002). Multifractal detrended fluctuation
% analysis of nonstationary time series. Physica A: Statistical Mechanics
% and its Applications, 316(1), 87-114.