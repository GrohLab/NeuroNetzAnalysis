classdef FourierSpectrum
% FOURIERSPECTRUM performs the 
    properties (SetAccess = 'private')
        Spectrum
        Fx
        Base2 = 2;
    end
    properties (Dependent)
        Magnitude
        Angle
        Ns
    end
    methods
        function obj = FourierSpectrum(cw)
            if nargin == 1 && isa(cw,'ContinuousWaveform')
                obj.Base2 = floor(log2(cw.NSamples));
                Lower2Pwr = 2^obj.Base2;
                win = blackmanharris(Lower2Pwr)';
                temp = cw.Data(1:Lower2Pwr);
                temp = win.*temp;
                obj.Spectrum = fftshift(fft(temp));
                obj.Fx = -cw.SamplingFreq/2:cw.SamplingFreq/Lower2Pwr:...
                    (cw.SamplingFreq*(Lower2Pwr - 2))/(2*Lower2Pwr);
            else
                % Blank initialization
                % Nothing to be done here...
                obj.Spectrum = [];
                obj.Fx = [];
            end
        end
        function magnitude = get.Magnitude(obj)
            magnitude = 20*log10(abs(obj.Spectrum));
        end
        function angl = get.Angle(obj)
            angl = unwrap(angle(obj.Spectrum));
        end
        function NumberOfSamples = get.Ns(obj)
            NumberOfSamples = exp(obj.Base2);
        end
        function p = plot(obj,varargin)
            figure;p(1) = subplot(2,1,1);plot(obj.Fx,obj.Magnitude,varargin{:});
            title('Magnitude');xlabel('Frequency [Hz]');ylabel('dB')
            p(2) = subplot(2,1,2);plot(obj.Fx,obj.Angle);
            title('Angle');xlabel('Frequency [Hz]');ylabel('Angle [rad]')
            linkaxes(p,'x');
        end
    end
end