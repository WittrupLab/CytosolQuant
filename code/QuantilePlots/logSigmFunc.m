function [y] = logSigmFunc(x, logEC50, HillSlope, Top, Bottom)
%logSigmFunc PRISM log(x) vs response sigmoidal
%   Detailed explanation goes here
y = Bottom + (Top-Bottom)./(1+10.^((logEC50-x).*HillSlope));
end

