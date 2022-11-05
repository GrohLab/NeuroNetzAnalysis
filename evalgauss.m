function prob = evalgauss(x,mu,sig)
%% EVALGAUSS Evaluate and create gaussian curve.
% evalgauss(X,MU,SIG) evaluates the data in column vector X for the
% gaussian curve explained by the given paramenters MU and SIG.
%
% $\aleph\left(x_i|\mu,\sigma\right)=
% \frac{1}{2\pi^{\frac{d}{2}}\sigma^\frac{1}{2}}
% e^{-\frac{\left(x_i-\mu\right)'\left(x_i-\mu\right)}{2\sigma}}$
%% Function
auxc = 1/sqrt(2*pi*sig);
auxe = ((x - mu).^2)/(2*sig);
prob = auxc * exp(-auxe);
end