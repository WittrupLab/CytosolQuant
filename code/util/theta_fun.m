function [f, f_int] = theta_fun(t, t0, tdelta, theta, scale)
%%
%   simple function for eval effect
%   
%   -exp(theta - min(0,t-t0-tdelta).^2/scale(1)^2) -max(0,t-t0 - tdelta ).^2/scale(2)^2 )
%
%
%%
scale = exp(scale);
tdelta = exp(tdelta);
f = -exp(theta - min(0,t-t0 - tdelta ).^2/scale(1)^2 - max(0,t-t0 - tdelta ).^2/scale(2)^2 ).* (t>=t0);


%%
%   trapz integration
%%
dt = t(2) - t(1);
f_int = f;
f_int(2:end-1) = 2*f_int(2:end-1);
f_int = dt * cumsum( f_int)/2;