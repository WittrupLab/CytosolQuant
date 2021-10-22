function [loss] = optim_theta(theta, t, A, y)
%%
%   simple function for eval effect
%   
%
%
%
%%
t0      = theta(1);
tdelta  = theta(2);
theta_f = theta(3);
scale   = theta(4:5);
[~, f_int] = theta_fun(t, t0, tdelta, theta_f, scale);
loss = sum((y - A * f_int).^2);