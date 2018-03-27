                            %OUTPUTS
function [ SE,...           Spectral Entropy
           ADRe,...         Entropy based alpha-delta ratio
           PRIe,...         Entropy based power ratio index
           Hj,...           Within band entropies
           beta,...         Global linear fit
           pfit,...         Distance to the perfect fit
           distLin ] =...   Goodness for the global and segmented fit
           spectEnt(...     %INPUTS
           mean_spectrum... Spectrum per channel
           )
%spectEnt Spectral entropy of the given signal.
%   Given a window size in ms with a certain overlap in percentage, the
%   function computes the power spectrum (PS) for the also given signal and
%   returns $\beta$ for the linear fit $log(PS)\alpha -\beta log(f)$.
%% Finding the paper described frequency band (0.5 - 35.0 Hz)
brain_waves = [ 7.0, 13.0;...   % Alpha
               13.1, 30.0;...   % Beta
                0.5,  4.0;...   % Delta
                4.1,  6.9;...   % Theta
                0.5, 30.0];     % Total
[idxs] =getIdxBrainW(mean_spectrum(1,:),brain_waves);
f_axis = mean_spectrum(1,idxs(5,1):idxs(5,2));

%% Computing the frequency complexity measures
[beta,pfit,distLin] = PSDslope(f_axis,mean_spectrum,idxs(5,1):idxs(5,2));
[SE,ADRe,PRIe,Hj] = spectral_Entropy(mean_spectrum,idxs);

end
                            %%OUTPUTS
function [beta,...          Global linear fitting
          feat,...          Distance to the perfect linear fit
          lin] =...         Goodness of global and segmented fit
          PSDslope(...      %%INPUTS
          freq_band,...     Frequency axis
          mean_spectrum,... Matrix with the signals' spectrums
          idx_f...          Index indicating 0.5 - 35.0 Hz range
          )
% PSDSLOPE Returns the linear fit for the log-power spectrum density vs.
% log-frequency. Power law.

%% Separating the slow and fast frequencies.
sb = freq_band <= 7;    %Frequency equal and below 7 Hz
fb = freq_band > 7;     %Frequency above 7 Hz
lf = log2(freq_band);   %Log-Frequency
lp = log2(mean_spectrum(2:end,idx_f));  %Log-Power
channels=size(lp,1);
%% Line Fitting
%Fitting three lines: a global, a low frequency one (delta and theta) and
%a fast frequency one (alpha and beta).
M = [lf' ones(length(lf),1)];
M1 = [lf(sb)' ones(length(lf(sb)),1)];
M2 = [lf(fb)' ones(length(lf(fb)),1)];
beta1 = pinv(M1) * lp(:,sb)';
beta2 = pinv(M2) * lp(:,fb)';
beta = pinv(M) * lp';
%Evaluating the polynomial of first order.
fx = beta(1,:)'*lf;
fx1 = beta1(1,:)'*lf;
fx2 = beta2(1,:)'*lf;
SS_tot = zeros(1,size(fx,1));
SS_tot1 = SS_tot;
SS_tot2 = SS_tot1;
%Computing the goodness of the fit.
for idx = 1:size(fx,1)
    fx(idx,:)=fx(idx,:)+beta(2,idx);
    fx1(idx,:)=fx1(idx,:)+beta1(2,idx);
    fx2(idx,:)=fx2(idx,:)+beta2(2,idx);    
    SS_tot(idx) = sumsqr(lp(idx,:)-mean(lp(idx,:)));
    SS_tot1(idx) = sumsqr(lp(idx,sb)-mean(lp(idx,sb)));
    SS_tot2(idx) = sumsqr(lp(idx,fb)-mean(lp(idx,fb)));
end
%% Goodness for the fit
%Computing the R squared for all the segments.
rsqr = 1 - (sum((lp - fx).^2,2)./SS_tot');
rsqr1 = 1 - (sum((lp(:,sb) - fx1(:,sb)).^2,2)./SS_tot1');
rsqr2 = 1 - (sum((lp(:,fb) - fx2(:,fb)).^2,2)./SS_tot2');
%Distance to the 'perfect fit'
unos = [1,1,1];             % R^2 in all segments--perfect fit.
lin = [rsqr,rsqr1,rsqr2];   % Real goodness
feat = distmatrix(unos,lin);% Distance to the perfect fit.
end
                                %%OUTPUTS
function [SE,...                Spectral entropy
          ADRe,...              Entropy based alpha-delta ratio
          PRIe,...              Entropy based power ratio index
          Hj] =...              Within band entropies
          spectral_Entropy(...  %%INPUTS
          mean_spectrum,...     Matrix with the signals' spectrums
          idx_f...              Index indicating 0.5 - 35.0 Hz range.
          )
% SPECTRAL_ENTROPY Returns the entropy for the given power spectrum
% density.

%% Power computing
PSD = 20*log10(mean_spectrum(2:end,idx_f(5,1):idx_f(5,2)));
for b=1:4
    Norms(b) = log2(idx_f(b,2)-idx_f(b,1)+1);
end
%% Entropy per channel
for ch = 1:size(PSD,1)
    rPSD(ch,:) = PSD(ch,:)./sum(PSD(ch,:));
    II(ch) = -rPSD(ch,:)*log2(rPSD(ch,:))';
%     Ha(ch) = -rPSD(ch,
    apwrs = rPSD(ch,(idx_f(1,1):idx_f(1,2))-idx_f(3,1)+1);
    bpwrs = rPSD(ch,(idx_f(2,1):idx_f(2,2))-idx_f(3,1)+1);
    dpwrs = rPSD(ch,(idx_f(3,1):idx_f(3,2))-idx_f(3,1)+1);
    tpwrs = rPSD(ch,(idx_f(4,1):idx_f(4,2))-idx_f(3,1)+1);
    % Computing the so called 'within-band' entropies
    Ha(ch) = - (apwrs/sum(apwrs))*log2(apwrs/sum(apwrs))';
    Hb(ch) = - (bpwrs/sum(bpwrs))*log2(bpwrs/sum(bpwrs))';
    Hd(ch) = - (dpwrs/sum(dpwrs))*log2(dpwrs/sum(dpwrs))';
    Ht(ch) = - (tpwrs/sum(tpwrs))*log2(tpwrs/sum(tpwrs))';
    % Computing the 'between-band' entropies
    Qa(ch) = - (sum(apwrs)*log2(sum(apwrs)));
    Qb(ch) = - (sum(bpwrs)*log2(sum(bpwrs)));
    Qd(ch) = - (sum(dpwrs)*log2(sum(dpwrs)));
    Qt(ch) = - (sum(tpwrs)*log2(sum(tpwrs)));
    Q(ch,:) = [sum(apwrs),sum(bpwrs),sum(dpwrs),sum(tpwrs)];
    Hj(ch,:) = [Ha(ch),Hb(ch),Hd(ch),Ht(ch)];
    for j = 1:4
        I(ch,j) = -(Q(ch,j)*log2(Q(ch,j))) + Q(ch,j)*Hj(ch,j);
        P(ch,j) = 100* (I(ch,j)/II(ch));
    end
end
for b=1:4
    Hj(:,b)=Hj(:,b)/Norms(b);
end
SE = II./log2(length(PSD));
ADRe = Hj(:,1)./Hj(:,3);
PRIe = (Hj(:,3)+Hj(:,4))./(Hj(:,1)+Hj(:,2));
end







