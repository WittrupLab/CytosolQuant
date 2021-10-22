function [eps_inj,sigma_inj, mcmc_obj] = sample_injection(eps_inj, mu, iSigmas, theta_inj, mcmc_obj, Y, inj)
%%
%   Samples inj^* given field parame
%
%   mu         - the mean of the field
%   iSigma     - precision matrix
%   theta_inj  - parameter for injection function
%   inj_param  - sigma
%   mcmc_obj   -  mcmc sampling
%
%%



ind = ~isnan(Y);
[~, T] = size(Y);


t0          = theta_inj(1);
delta       = exp(theta_inj(2));
beta        = exp(theta_inj(3:4));
scale       = exp(theta_inj(5:6));

dt=1/10;


alpha_gamma = T/2+1;
beta_gamma = sum(eps_inj.^2)/2 + 0.001;
sigma_inj = 1;%sqrt(1./gamrnd(alpha_gamma,1/beta_gamma));
mu_inj = sigma_inj^2;
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
  [f, f_int] = injection_func(time, inj(i)*exp(eps_inj(i)), t0, delta, beta, scale, dt);
  Ytb = Y(ind(:,i),i) - mu - A * f_int;
  
  lik_old = - eps_inj(i) - 0.5 * (eps_inj(i)- mu_inj)^2/sigma_inj^2;
  lik_old = lik_old - 0.5*(Ytb'*iSigmas{i}*Ytb);
  
  eps_inj_star =  eps_inj(i) + mcmc_obj.sigma(i) * randn;
  
  lik_star = - eps_inj_star - 0.5 * (eps_inj_star - mu_inj)^2/sigma_inj^2;
  [f, f_int] = injection_func(time, inj(i)*exp(eps_inj_star), t0, delta, beta, scale, dt);
  Ytb = Y(ind(:,i),i) - mu - A * f_int;
  lik_star = lik_star - 0.5*(Ytb'*iSigmas{i}*Ytb);
  
  if log(rand) < lik_star - lik_old
     eps_inj(i)  = eps_inj_star;
     mcmc_obj.acc_sum(i)         = mcmc_obj.acc_sum(i) + 1;
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
