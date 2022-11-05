function num = bin2num(binArray, num)
%BIN2NUM converts a binary number from left to right.
%   The function rises two to the power of the length of the binArray - 1
%   if it contains true and deletes the value from the array for the next
%   computation.
%           num = bin2num(binArray, num)
%       INPUTS
%           binArray - vector containing the logical number
%       OUTPUT
%           num - the cumulative variable for the final conversion
Nb = length(binArray);
if Nb > 0
    num = num + (2^(Nb - 1)*binArray(1));
    binArray(1) = [];
    num = bin2num(binArray, num);
end
end