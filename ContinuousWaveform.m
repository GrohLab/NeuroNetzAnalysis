classdef ContinuousWaveform < GeneralWaveform
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Spectrum FourierSpectrum;
    end
    properties (Dependent)
        TMean
    end
    
    methods
        function obj = ContinuousWaveform(inputArg1,inputArg2)
            %CONTINUOUSWAVEFORM creates an instance of the
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function trigMean = get.TMean(obj)
            %METHOD1 This method should take into consideration the
            %timestamps of the DiscreteWaveform class and get an average of
            %the continuous waveform.
            %   Detailed explanation goes here
            trigMean = obj.Data;
        end
         function spectrum = get.Spectrum(obj)
            spectrum = FourierSpectrum(obj);
        end
    end
end

