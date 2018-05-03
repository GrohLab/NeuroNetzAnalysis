classdef BurstTrain
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = 'private')
        First = [];
        Count = [];
        Length = [];
    end
    
    methods
        function obj = BurstTrain(spObj,isi)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 0
                obj.First = 1;
                obj.Count = 0;
                obj.Length = 0;
            elseif nargin == 1
                obj.First = spObj.Spikes;
                obj.Count = zeros(spObj.Ns);
                obj.Length = obj.Count;
            elseif nargin == 2
                
            else
                
            end
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
    methods (Static)
        function isiEst = isiDesicion(spObj,verbose)
            lisi = log(diff(spObj.Spikes)/spObj.samplingFreq);
            if nargin == 1
                verbose = true;
            end
            prms = emforgmm(lisi,5,1e-7,0);
            thrshs = findthreshGMM(prms,lisi);
            isiEst = thrshs(1);
            if verbose
                figure;histogram(lisi)
            end
            
        end
    end
end

