function [signal_out] = iir50CombFilter( signal, sampling_frequency, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;

addRequired(p, 'signal', @(x) isvector( x ) && isnumeric( x ) );
addRequired(p, 'sampling_frequency', @(x) isnumeric( x ) && x > 0 )
addParameter(p, 'Q', 35, @(x) isnumeric( x ) && isscalar( x ) && x > 0 );

parse(p, signal, sampling_frequency, varargin{:} );

signal = p.Results.signal;
sampling_frequency = p.Results.sampling_frequency;
Q = p.Results.Q;

fnOpts = {'UniformOutput', false};

wo = 50:50:(sampling_frequency/2);
[b, a] = arrayfun(@(w) iirnotch( w, (2*w)/Q ... Notch bandwidth using a 'Q' factor.
    ), (2*wo)/sampling_frequency, fnOpts{:} );

fprintf( 1, 'Filtering the signal for 50 Hz and its harmonics...' )
signal_out = signal;
for cfb = 1:numel(b)
    signal_out = filtfilt( b{cfb}, a{cfb}, signal_out );
end
fprintf(' done!\n');

end