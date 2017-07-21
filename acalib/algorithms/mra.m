function [approx] = mra(data, wavelet, levels)
%It performs mra through the 3D DWT
n = levels;

%Compute the wavelet decomposition of the 3D data
WT = wavedec3(data,levels,wavelet);

%Reconstruct from coefficients the approximations for each level
approx = cell(1,levels);
for k = 1:levels
    % Compute reconstructed approximation, i.e. the low-pass component.
    approx{k} = waverec3(WT,'a',k);  
end
clear('WT');
end

