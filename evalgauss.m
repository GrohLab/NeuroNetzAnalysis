function prob = evalgauss(x,mu,sig)
%% EVALGAUSS Evaluate and create gaussian curve.
% evalgauss(X,MU,SIG) evaluates the data in column vector X for the
% gaussian curve explained by the given paramenters MU and SIG.
%
% $\aleph\left(x_i|\mu,\Sigma\right)=
% \frac{1}{2\pi^{\frac{d}{2}}\Sigma^\frac{1}{2}}
% e^{-\frac{\left(x_i-\mu\right)'\left(x_i-\mu\right)}{2\Sigma}}$
%% Function
d=1;
auxc = 1/(sqrt((2*pi)^d)*sqrt(sig));
auxe = ((x - mu).^2)/(2*sig);
prob = auxc * exp(-auxe);
end