classdef gmmpdf
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = 'protected', SetAccess = 'private')
        parameters = [];
    end
    properties (GetAccess = 'private', SetAccess = 'protected')
        originalSum = [];
    end
    properties
        dataDomain = [];
        resolution = [];
    end
    
    properties (Dependent)
        pdf
    end
    
    methods
        function obj = gmmpdf(varargin)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
           for narg = 1:nargin
               cols = size(varargin{narg},2);
               if cols == 3
                   obj.parameters = varargin{narg};
               elseif cols == 2
                   obj.dataDomain = varargin{narg};
               elseif cols == 1
                   % Must validate the resolution or bining of the
                   % y-values.
                   obj.resolution = varargin{narg};
               else
                   warning('Input %d is not recognized. Ignoring...\n',narg) 
               end
           end
           if isempty(obj.dataDomain)
               if ~isempty(obj.parameters)
                   % disp('Estimating the pdf domain given its parameters:')
                   [~,GIdx] = max(obj.parameters(:,2));
                   [~,sIdx] = min(obj.parameters(:,2));
                   bigSTD = max(obj.parameters(:,3));
                   x_low = obj.parameters(sIdx,2) - bigSTD*5;
                   x_high = obj.parameters(GIdx,2) + bigSTD*5;
                   obj.dataDomain = [x_low, x_high];
                   % fprintf('X: (%f, %f)\n',x_low,x_high)
               else
                   error(['A probability density function for a GMM ',...
                       'cannot be computed without parametrization! ',...
                       'Aborting...\n'])
               end
           end
           if isempty(obj.resolution)
               obj.resolution = diff(obj.dataDomain)/((2^10)-1);
           end
        end
        
        function pdf = get.pdf(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            % Section thought for 'un-normalization'
            %[pdf, OSum] = genP_x(obj.parameters,obj.getDataDomain());
            %if isempty(obj.originalSum)
            %    obj.originalSum = OSum;
            %end
            pdf = genP_x(obj.parameters,obj.getDataDomain());
        end
        
        function H = getEntropy(obj)
            obj.pdf = obj.genPDF();
            N = length(obj.pdf);
            if sum(obj.pdf)~=1
                obj.pdf = obj.pdf/sum(obj.pdf);
            end
            H = -dot(obj.pdf,log(obj.pdf))/log(N);
        end
        
        function dataAxis = getDataDomain(obj)
            dataAxis = obj.dataDomain(1):obj.resolution:obj.dataDomain(2);
        end
        
        function plotPDF(obj)
            figure('name','Probability density function');
            plot(obj.getDataDomain(),obj.genPDF());grid('on')
        end
        
        function kld = comparePDF(obj,obj2)
            if isa(obj2,'gmmpdf')
                if isempty(obj.pdf)
                    obj.genPDF();
                end
                if isempt(obj2.pdf)
                    obj2.genPDF();
                end
                if length(obj.pdf) ~= length(obj2.pdf)
                    obj2.dataDomain = obj.dataDomain;
                    obj2.resolution = obj.resolution;
                    obj2.genPDF();
                end
                kld = dot(obj.pdf,log(obj.pdf ./ obj2.pdf))/log(...
                    length(obj.pdf));
            end
        end
    end
end

