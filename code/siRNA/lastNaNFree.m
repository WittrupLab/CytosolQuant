function [index] = lastNaNFree(array)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
a=1:length(array);
index=min(a(isnan(array)))-1;
end



