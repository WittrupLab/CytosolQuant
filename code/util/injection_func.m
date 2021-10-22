function [f,f_int] = injection_func(t, inj, t0, delta, beta, scale, dt)
%%
%   The injection function 
%   f(t) = - beta_1 * inj^{\beta_2}e^{ -(t-t0-delta)^2/scale} I(t>=t0);
%
%
%
%%
if nargin < 7
    dt = 1;
end
index = t>=t0;
t_mod = t(index)- t0 - delta;
f = zeros(size(t));
f(index) = - beta(1) * inj^(beta(2)) .*  exp( -min(0,t_mod).^2/scale(1)^2 - max(0,t_mod).^2/scale(2)^2);

%%
%   trapz integration
%%
f_int = f;
f_int(2:end-1) = 2*f_int(2:end-1);
f_int = dt * cumsum( f_int)/2;