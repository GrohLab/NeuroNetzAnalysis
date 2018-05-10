classdef FourierSpectrum
    properties (SetAccess = 'private')
        Spectrum
        Fx
    end
    properties (Dependent)
        Magnitude
        Angle
    end
    methods
        function obj = FourierSpectrum(gw)
            exptwo = floor(log2(gw.NSamples));
            desSize = 2^exptwo;
            win = blackmanharris(desSize)';
            temp = gw.Data(1:desSize);
            temp = win.*temp;
            obj.Spectrum = fftshift(fft(temp));
            obj.Fx = -gw.SamplingFreq/2:gw.SamplingFreq/desSize:...
                (gw.SamplingFreq*(desSize - 2))/(2*desSize);
        end
        function magnitude = get.Magnitude(obj)
            magnitude = 20*log10(abs(obj.Spectrum));
        end
        function angl = get.Angle(obj)
            angl = unwrap(angle(obj.Spectrum));
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