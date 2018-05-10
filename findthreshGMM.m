function thrsh = findthreshGMM(params,data)
%FINDTHRESHGMM finds the peaks and vallies of the probability density
%function given its estimated parameters and the data from which the
%parameters were estimated.
data = data(:);
xaxis = linspace(min(data),max(data),numel(data));
p_x = genP_x(params,xaxis); % Final distributions
% The derivative of the logarithm of the absolute values from a first
% derivative returns spikes where the first derivative is equal to zero.
ad_px = diff(log(abs(diff(log(p_x))))); 
posPts = ad_px > 0.5;
negPts = ad_px < -0.5;
% Need to take care when the positive and the negative spikes are not equal
thrsh = (xaxis(diff(posPts)==-1) + xaxis(diff(negPts)==1))/2;
% DELETED SECTION 2
end


%% DELETED SECTION 2
% near_zero = ((max(ad_px) - min(ad_px))*1/(2^cero)) + min(ad_px);
% zcs = ad_px < near_zero;            % Finding segments which are near zero
% bgn = xaxis(diff(zcs)==1);          % Beginning of the ones segment
% fin = xaxis(diff(zcs)==-1);         % End of the ones segment
% pav = min(length(bgn),length(fin)); % Peaks And Vallies
% zcs = (bgn(1:pav)+fin(1:pav))/2;    % Middle point of the segment
% if ~isempty(zcs)
%     thrsh = zcs;
% else
%     thrsh = NaN;
% end