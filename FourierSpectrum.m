classdef FourierSpectrum < handle
% FOURIERSPECTRUM performs the 
    properties (SetAccess = 'private')
        FourierTransform
        Fx
        Magnitude
        Phase
        N
        PowerSpectrum
        Fs
    end
    
    methods        
        % Class constructor
        function obj = FourierSpectrum(input_signal, fs)
            if isnumeric(input_signal)
                % The signal will be zero padded and windowed using a
                % _specific_window_. A question arises: is it necessary to
                % compute a sonogram or only with the direct Fourier
                % transform is enough?
                [Nch, Ns] = size(input_signal);
                obj.Fs = fs;
                if Nch > Ns
                    input_signal = input_signal';
                    Ns2 = Nch;
                    Nch = Ns;
                    Ns = Ns2;
                end
                % Normalizing the amplitude of the signal (the maxima
                % should be equal to 1 or -1)
                % Removal of the DC  (0 Hz) component
                input_signal = input_signal - mean(input_signal);
                % Amplitude normalization
                input_signal = input_signal / rms(input_signal);
                % Determination for a sensible padding.
                N1 = floor(log2(Ns));
                N2 = nextpow2(Ns);
                if N1 ~= N2
                    d1 = abs((2^N1) - Ns)/(2^N2 - 2^N1);
                    d2 = abs((2^N2) - Ns)/(2^N2 - 2^N1);
                    if d1 > d2
                        % Zero padding
                        obj.N = 2.^N2;
                        padSz = (obj.N-Ns);
                        if mod(padSz,2)
                            pad_signal = padarray(input_signal,[0,padSz],'post');
                        else
                            pad_signal = padarray(input_signal,[0,padSz/2]);
                        end
                    else
                        % Cutting the signal to the immediate lower power of 2
                        obj.N = 2.^N1;
                        padSz = Ns - obj.N;
                        if ~mod(padSz,2)
                            pad_signal = input_signal(padSz/2 + 1:...
                                end - padSz/2);
                        else
                            pad_signal = input_signal(1:end - padSz);
                        end
                    end
                elseif (N1 - log2(Ns)) == 0
                    % No padding needed
                    pad_signal = input_signal;
                    obj.N = Ns;
                end
                % Window the padded signal
                auxWin = hann(obj.N)';
                try
                    win_signal = pad_signal .* auxWin;
                catch
                    disp('Not so much memory left! Using single samples')
                    win_signal = zeros(size(pad_signal),'single');
                    for cs = 1:obj.N-4
                        % Vectorization of independent samples
                        win_signal(cs) = pad_signal(cs) * auxWin(cs);
                        win_signal(cs+1) = pad_signal(cs+1) * auxWin(cs+1);
                        win_signal(cs+2) = pad_signal(cs+2) * auxWin(cs+2);
                        win_signal(cs+3) = pad_signal(cs+3) * auxWin(cs+3);
                    end
                end
                % Compute the Fourier Transformation of the _input_signal_
                % together with the frequency axis. The factor of 2
                % compensates the Hann windowing.
                obj.FourierTransform = 2 * fftshift(fft(win_signal));
                obj.Fx = -fs/2:fs/obj.N:fs/2 - fs/obj.N;
            else
                disp('Not implemented yet')
            end
        end
        
        %% GET functions
        function PSD = get.PowerSpectrum(obj)
            % The power spectrum is computed according to the energy for
            % each bin size $\frac{Fs}{N}$ 
            %PSD = (2*obj.FourierTransform.*conj(obj.FourierTransform))...
            %    /(obj.N^2);
            PSD = (1/(obj.Fs*obj.N))*abs(obj.FourierTransform).^2;
        end
        
        function ph = get.Phase(obj)
            ph = angle(obj.FourierTransform);
        end
        
        function mg = get.Magnitude(obj)
            mg = abs(obj.FourierTransform);
        end
        
        %% Half signals
        function posFourier = getHalFourier(obj)
            posFourier = FourierSpectrum.getRightHalf(...
                obj.FourierTransform);
            posFourier(2:end) = 2 * posFourier(2:end);
        end
        
        function posPower = getHalfPower(obj)
            posPower = FourierSpectrum.getRightHalf(...
                obj.PowerSpectrum);
            posPower(2:end) = 2 * posPower(2:end);
        end
        
        function halPh = getHalPhase(obj)
            halPh = FourierSpectrum.getRightHalf(obj.Phase);
        end
        
        function halFx = getHalFx(obj)
            halFx = FourierSpectrum.getRightHalf(obj.Fx);
        end
        %% Normalization functions (obsolete)
        function outSignal = normalizeSum(obj)
            outSignal = obj.PowerSpectrum/sum(obj.PowerSpectrum);
        end
        
        function outSignal = normalizeLength(obj)
            outSignal = obj.PowerSpectrum/obj.N;
        end
        
        function outSignal = normalizeDeltaF(obj)
            outSignal = obj.PowerSpectrum/(obj.Fs/obj.N);
        end
        %% Plotting functions
        function p = plotFrequencyDomain(obj,varargin)
            p = plot(obj.Fx,obj.PowerSpectrum,varargin{:});
            ylabel('Power [\muV^2]');xlabel('Frequency [Hz]')
            title('Power Spectrum')
        end
        
        function p = plotFrequency_dB(obj,varargin)
            p = plot(obj.Fx,10*log10(obj.PowerSpectrum),varargin{:});
            ylabel('dB');xlabel('Frequency [Hz]')
            title('Power Spectrum')
        end
        
        function p = plotHalfFrequencyDomain(obj,varargin)
            p = plot(obj.getHalFx,obj.getHalfPower,varargin{:});
            ylabel('Power [\muV^2]');xlabel('Frequency [Hz]')
            title('Power Spectrum')
        end
        
        function p = plotHalfFrequency_dB(obj,varargin)
            p = plot(obj.getHalFx,10*log10(obj.getHalfPower),varargin{:});
            ylabel('dB');xlabel('Frequency [Hz]')
            title('Power Spectrum')
        end
    end
    methods (Static, Access = 'private')
        function halfSignal = getRightHalf(wholeSignal)
            N = length(wholeSignal);
            if mod(N,2)
                halfSignal = wholeSignal(ceil(N/2):end);
            else
                halfSignal = wholeSignal((N/2)+1:end);
            end
        end
    end
end