function [index] = maxNonNaN(array)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
a=1:length(array);
index=max(a(not(isnan(array))));
end

