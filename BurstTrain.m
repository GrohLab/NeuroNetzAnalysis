classdef BurstTrain
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = 'private')
        First = [];
        Count = [];
        Length = [];
        ISI = 1e-3;    % 1 ms. The units are always seconds otherwise indicated
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
                obj.ISI = BurstTrain.isiDesicion(spObj,0);
                
                obj.First = spObj.Spikes;
                obj.Count = zeros(spObj.Ns);
                obj.Length = obj.Count;
            elseif nargin == 2
                obj.ISI = isi;
                
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
            lisi = log(diff(spObj.Spikes)/spObj.SamplingFreq);
            if nargin == 1
                verbose = true;
            end
            pisi = emforgmm(lisi,8,1e-6,1);
            thrsh = findthreshGMM(pisi,lisi);
            pdfVals = genP_x(pisi,thrsh);
            [~,thIdx] = min(pdfVals);
            isiEst = exp(thrsh(thIdx));
            if verbose
                xaxis = linspace(min(lisi),max(lisi),numel(lisi));
                xaxisT = exp(xaxis);
                figure;semilogx(xaxisT,genP_x(pisi,xaxis))
            end
        end
        function [bIdx, tIdx, spIdx] = getInitialBurstSpike(spObj,maxISI)
            % GETINITIALBUSRTSPIKE gets the first spike time for each burst. The tonic
            % spikes are unmodified according to the maximum inter-spiking-interval.
%             if ~isrow(spObj.Spikes)
%                 spObj.Spikes = spObj.Spikes';
%             end
            delta_spT = cat(find(size(spObj.Spikes)~=1),inf,diff(spObj.Spikes));
            % delta_spT = [inf,diff(spT)];
            fsIdx = delta_spT >= maxISI;    % First spike index (burst-wise)
            bsIdx = [~fsIdx(2:end) & fsIdx(1:end-1),fsIdx(end)];
            tsIdx = ~bsIdx & fsIdx;
            bIdx = spObj.Spikes(bsIdx);
            tIdx = spObj.Spikes(tsIdx);
            spIdx = spObj.Spikes(bsIdx | tsIdx);
        end
    end
end

