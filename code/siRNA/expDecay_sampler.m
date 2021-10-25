function [theta, acc, obj] = expDecay_sampler(theta, i, obj)
%%
%
%
%   
%   obj - 
%           .Sigma
%           .theta_hat
%           .mcmc_type
%           .sigma_mcmc
%%
[R,p] = chol(obj.Sigma);
if p > 0
R = diag(sqrt(diag(obj.Sigma)));    
end

if obj.mcmc_type == 1
    sigma = min(1, obj.sigma_mcmc);
    w_theta = (1- sigma^2)^(1/2);
    theta_mu   = w_theta * obj.theta_hat      + (1- w_theta) * theta;
    Z_s = randn(size(theta));
    theta_star = theta_mu + obj.sigma_mcmc  * R' * Z_s;
    theta_mus  = w_theta * obj.theta_hat + (1-w_theta) * theta_star;
else
    theta_mu  = theta;
    Z_s = randn(size(theta));
    theta_star = theta_mu + obj.sigma_mcmc  * R' * Z_s;
    theta_mus = theta_star; 
end

Zs = R'\( theta_star - theta_mu);
lq_xs = - 1/(2 * obj.sigma_mcmc^2) * (Zs' * Zs);
Z0 = R'\( theta - theta_mus);
lq_x0 = - 1/(2 * obj.sigma_mcmc^2) * (Z0' * Z0);

acc_p = obj.logpi(theta_star) - obj.logpi(theta) + lq_x0 - lq_xs;
acc = 0;
if log(rand) < acc_p
   acc =1;
   theta = theta_star;
end
% n?r man b?rjar anv?nda Sigma ?ndra sigma till
% sigma/min(diag(sqrt(Sigma)))
w = 1/i;

obj.theta_hat = (1/w-1) * w * obj.theta_hat + w * theta;
obj.Sigma_hat = (1/w-1) * w * obj.Sigma_hat + w * (theta - obj.theta_hat) * (theta - obj.theta_hat)';
