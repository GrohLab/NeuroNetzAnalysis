function Pstruct = getPDFFromGMM(params,xaxis)
%GETPDFFROMGMM returns a structure with the relevant information about the
%associated PDF for the given data and the parameters from a Gaussian
%Mixture Model (GMM).
%   Detailed explanation goes here

% Redifine the PDF structure. Use use defined X axis. Same domain for every
% PDF.

pdf = genP_x(params,xaxis);
domain = [min(axis), max(axis)];
Pstruct = struct('PDF',pdf,'Axis',axis,'Domain',struct('Min',domain(1),...
    'Max',domain(2)),'Entropy',getEntropyFromPDF(pdf));
Pstruct.Integrate =...
    @(x_low,x_high)...
    (sum(Pstruct.PDF(find(Pstruct.Axis >= x_low,1):...
    find(Pstruct.Axis >= x_high)-1),'omitnan'));
Pstruct.Normalize = @()(Pstruct.PDF/sum(Pstruct.PDF));
Pstruct.cutAndDownsample = @(x_low,x_high,N)...
    (Pstruct.PDF(round(linspace(find(Pstruct.Axis >= x_low,1),...
    find(Pstruct.Axis >= x_high,1)-1,N))));
Pstruct.CDF = @() (cumsum(Pstruct.PDF));
end