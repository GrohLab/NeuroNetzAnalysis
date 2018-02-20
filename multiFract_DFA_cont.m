                                %OUTPUTS
function [Ht,...                Time depending Hurst exponents
          Hb,...                X-axis for the pdf
          Pdf,...               Hurst probability distribution
          mfs]...               Multifractal Spectrum
    = multiFract_DFA_cont(...   %INPUTS
          signal,...            Time-series to be analyzed
          ts,...                Time axis
          scale_b,...           Monofractal scales
          scale_s,...           Multifractal scales
          m,...                 m-order polynomial to detrend
          econ)...              Economic estimation of the PDF
    %MULTIFRACT_DFA Computation of the multifractal spectrum using $q$-orders
%   Detailed explanation goes here

% Determining the profile as in Kantelhardt et al. 2002 (Integrating
% without the mean).

display('Computing monofractal analysis')
%% Monofractal analysis with q = 0
[profile,~,~,~,~,Fq0] = multiFract_DFA(signal,ts,scale_b,0,m);

%% Computing the fluctuation functions for every moving window
hmx = floor(scale_s(end)/2);
timeIdx = hmx+1:length(profile)-hmx;
for s = 1:length(scale_s)
    hseg = floor(scale_s(s)/2);
    for v = hmx+1:length(profile)-hmx
        % Computing the variation of the signal for each segment size using
        % a moving window centerd in v.
        idxs = v-hseg:v+hseg;
        [polft,~]=detrend_profile(m,ts(idxs),profile(idxs));
        bs{s}(v) = sqrt(mean((profile(idxs)-polft).^2));
    end
    % Avoiding the infinity
    F(s)=exp(0.5*mean(log(bs{s}(bs{s}~=0).^2)));
end
% Computing the hurst exponent from the monofractal fluctuation function
[~,ft] = detrend_profile(1,log2(scale_b),log2(Fq0));
H = ft(1);
avfit = log2(scale_s).*H + ft(2);
N = length(timeIdx);
%% Computing the Hurst exponents for each calculated variation
for s = 1:length(scale_s)
    bsC = bs{s}(timeIdx);
    if sum(bsC == 0)
        aux = sort(bsC(bsC~=0));
        bsC(bsC==0) = aux(1);
    end
    res_bs=avfit(s)-log2(bsC);
    lg = log2(N)-log2(scale_s(s));
    Ht(s,:)=res_bs./lg + H;
end

%% Computing the Multifractal Spectrum
Hrow = Ht(:);
aux=-log(mean(scale_s));
if econ
    nBins = ceil(length(Hrow)^(1/sqrt(2)));     % Estimating bins
    [f_,Hb]=hist(Hrow,nBins);                   % Computing the frequencies
    Pdf = f_/sum(f_);                           % Forcing a pdf
else
    [P,~]=emforgmm(Hrow,8,1e-3,0);              % Estimating a GMM
    [Pdf,Hb]=genP_x(P,Hrow);                    % Computing the PDF
end
Pdf_n = Pdf/max(Pdf);           % Normalizing for the logarithm
mfs = 1 - (log(Pdf_n)./aux);    % Multifractal spectrum
end