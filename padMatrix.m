function [padMat] = padMatrix(matx, dim, verb)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Verbose
if ~exist('verb', 'var')
    verb = false;
end
% Dimension verification
if ~(dim == 1 || dim == 2)
    fprintf(1, 'Error: %d dimension not supported!\n', dim)
    fprintf(1, 'Setting dimension to 1\n')
    dim = 1;
end
% User feedback
[Nr, Nc] = size(matx);
if verb && (dim == 1 || dim == 2)
    fprintf(1, 'Padding along the %d dimension\n', dim)
end
% Dimension workout
if dim == 1
    N = Nr;
else
    N = Nc;
end

bohWin = reshape(bohmanwin(N),([-dim+2, 0; 0, dim-1]*([N;N]-1))'+1);
matxWin = matx .* bohWin;

padDiff = (2^nextpow2(N)) - N;
padN = padDiff/2;
addOne = false;
if mod(padDiff, 2)
    padN = (padDiff - 1)/2;
    addOne = true;
end
pad1 = padarray(matxWin, [-dim+2, 0; 0, dim-1]*[padN;padN],0, 'pre');
padMat = padarray(pad1, [-dim+2, 0; 0, dim-1]*([padN;padN] + addOne), 0, 'post');

end

