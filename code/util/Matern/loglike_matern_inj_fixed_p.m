function lik = loglike_matern_inj_fixed_p(theta, iSigmas, theta_inj, Y, inj)
%%
%   Assumes that the parameters for the Sigma is fixed
%   theta - (5 x 1)
%           mu      - mean
%           kappa   - range
%           sigma   - standard dev
%           nu      - differntiablility
%           sigma_e - measurement error
%
%%



ind = ~isnan(Y);
[~, T] = size(Y);

mu      = theta(1);

t0          = theta_inj(1);
delta       = exp(theta_inj(2));
beta        = exp(theta_inj(3:4));
scale       = exp(theta_inj(5:6));

dt=1/10;

lik = 0;


for i=1:T
    
  f_index    = find(ind(:,i));
  y_start  = min(f_index);
  y_end    = max(f_index);
  N = length(y_start:y_end);
  N_x = ceil(N/dt);
  Yobs =  1 + floor((0:(N-1))/dt);
  A = speye(N_x);
  A = A(Yobs,:);
  A = A(ind(y_start:y_end,i),:);  
  time     = (y_start + (0:(N_x-1))*dt)';
  [f, f_int] = injection_func(time, inj(i), t0, delta, beta, scale, dt);
  Ytb = Y(ind(:,i),i) - mu - A * f_int;
  
  lik = lik - 0.5*(Ytb'*iSigmas{i}*Ytb);
  
%     if(i==10)
%       
%        subplot(211)
%        plot(time, f)
%        subplot(212)
%        plot(time, f_int)
%        drawnow
%        
%     end
end
