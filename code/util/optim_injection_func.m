function [loss] = optim_injection_func(theta, t, A, y, siRNA, logscale)
%%
%   simple function for eval effect
%   
%   theta - parameters
%   t     - observations scale
%   A     - observation matrix
%   y     - (N x k) 
%   siRNA - k x 1
%%
dt = t(2) - t(1);
t0      = theta(1);
tdelta  = exp(theta(2));
beta    = exp(theta(3:4));
scale   = exp(theta(5:6)); 
loss = 0;
for i=1:size(y,2)
    [~,f_int] = injection_func(t, ...
                               siRNA(i),...
                               t0, ...
                               tdelta, ...
                               beta, ...
                               scale,...
                               dt);
                           
  if(logscale)
    loss =loss +  nansum((y(:,i) - A * f_int).^2);
  else
    loss =loss +  nansum((y(:,i) - exp(A * f_int)).^2);
  end
end