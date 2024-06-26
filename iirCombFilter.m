function [signal_out] = iirCombFilter( signal, sampling_frequency, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;

addRequired(p, 'signal', @(x) isvector( x ) && isnumeric( x ) );
addRequired(p, 'sampling_frequency', @(x) isnumeric( x ) && x > 0 )
addParameter(p, 'Q', 35, @(x) isnumeric( x ) && isscalar( x ) && x > 0 );
addParameter(p, 'W0', 50, @(x) isnumeric( x ) && isscalar( x ) && ...
    x > 0 && x < (sampling_frequency/2) );
addParameter(p, 'verbose', true, @(x) islogical( x ) && isscalar( x ) );

parse(p, signal, sampling_frequency, varargin{:} );

signal = p.Results.signal;
sampling_frequency = p.Results.sampling_frequency;
Q = p.Results.Q;
w0 = p.Results.W0;
verbose = p.Results.verbose;

fnOpts = {'UniformOutput', false};

N_notch = sampling_frequency / (2*w0);

wo = (1:floor( N_notch) )*w0;
[b, a] = arrayfun(@(w) iirnotch( w, (2*w)/Q ... Notch bandwidth using a 'Q' factor.
    ), (2*wo)/sampling_frequency, fnOpts{:} );

if verbose
    fprintf( 1, 'Filtering the signal for 50 Hz and its harmonics...' )
end

signal_out = signal;
for cfb = 1:numel(b)
    signal_out = filtfilt( b{cfb}, a{cfb}, signal_out );
end

if verbose
    fprintf(' done!\n');
end

end