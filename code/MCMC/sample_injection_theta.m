function [eps_inj, mcmc_obj] = sample_injection_theta(sim, mu, iSigmas, theta_inj, mcmc_obj, Y, inj)
%%
%   Samples inj^* given field parame
%
%  eps_inj
%
%%



ind = ~isnan(Y);


t0          = theta_inj(1);
delta       = exp(theta_inj(2));
beta        = exp(theta_inj(3:4));
scale       = exp(theta_inj(5:6));

dt=1/10;

eps_sum_vec = zeros
sigma_inj = 1;%sqrt(1./gamrnd(alpha_gamma,1/beta_gamma));
mu_inj = sigma_inj^2;
for i=1:sim
    
  f_index    = find(ind);
  y_start  = min(f_index);
  y_end    = max(f_index);
  N = length(y_start:y_end);
  N_x = ceil(N/dt);
  Yobs =  1 + floor((0:(N-1))/dt);
  A = speye(N_x);
  A = A(Yobs,:);
  A = A(ind(y_start:y_end,i),:);  
  time     = (y_start + (0:(N_x-1))*dt)';
  [f, f_int] = injection_func(time, exp(eps_sum(1)), t0, delta, beta, scale, dt);
  Ytb = Y(ind(:,i),i) - mu - A * f_int;
  eps_inj = eps_sum(2);
  lik_old = - eps_inj(i) - 0.5 * (eps_inj- mu_inj)^2/sigma_inj^2;
  lik_old = lik_old - 0.5*(Ytb'*iSigmas{i}*Ytb);
  
  eps_inj_star =  eps_inj(i) + mcmc_obj.sigma(i) * randn;
  
  lik_star = - eps_inj_star - 0.5 * (eps_inj_star - mu_inj)^2/sigma_inj^2;
  [f, f_int] = injection_func(time, inj(i)*exp(inj_star + eps_inj_star), t0, delta, beta, scale, dt);
  Ytb = Y(ind(:,i),i) - mu - A * f_int;
  lik_star = lik_star - 0.5*(Ytb'*iSigmas{i}*Ytb);
  
  if log(rand) < lik_star - lik_old
     eps_sum_vec(i)  = eps_inj_star;
     acc_sum     = mcmc_obj.acc_sum(i) + 1;
  end
  
  [mcmc_obj.sigma(i), mcmc_obj.acc_sum(i)]  = AMCMC_RR(mcmc_obj.sigma(i), ...
                                                       mcmc_obj.acc_sum(i), ...
                                                       mcmc_obj.iter, ...
                                                       50, ...
                                                       0.4);
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
