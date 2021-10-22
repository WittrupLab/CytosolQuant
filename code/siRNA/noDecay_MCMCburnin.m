function theta = noDecay_MCMCburnin(ydata, xdata, sim)
%%
%   estimation parameters in the no expontial decay function
%  ydata - (N x 1) observations
%  xdata - (N x 1) location of the observations
%  sim   - (int)   number of simulations in the MCMC algorithm in the pre sampler 
%%
xdata = xdata(:);
ydata = ydata(:);

y_isnan = isnan(ydata);
xdata = xdata(y_isnan==0);
ydata = ydata(y_isnan==0);
burnin = round(0.5*sim);
obj = struct;
obj.Sigma_hat = zeros(3,3);
obj.theta_hat = zeros(3,1);
obj.Sigma     = eye(3); 
obj.sigma_mcmc = 0.001;
obj.mcmc_type = 2;

[argy, argm] = max(ydata(1:min(12,length(ydata))));
theta = [log(max(argm-0.5 - 1,1)),0,0]';
sigma = sqrt(var((ydata - noDecay(xdata,theta))));
acc_sum = 0;
for i= 1:sim
  obj.logpi = @(theta) - sum((ydata - noDecay(xdata,theta)).^2)/sigma^2 ;%- sum(theta.^2)/5^2;
 [theta, acc_i, obj] = expDecay_sampler(theta, i, obj);
 
 alpha = length(ydata)/2 + 1;
 beta   = sum((ydata - noDecay(xdata,theta)).^2)/2;
 sigma  = sqrt(1./gamrnd(alpha,1/beta));
 
    if acc_i==1
       acc_sum = acc_sum +1;
    end
    [obj.sigma_mcmc, acc_sum]  = AMCMC_RR(obj.sigma_mcmc, acc_sum, i, 50, 0.15);
    
    if i == burnin
        obj.sigma_mcmc = sqrt( min(diag(obj.Sigma)) *  obj.sigma_mcmc^2./min(diag(obj.Sigma_hat)));
        if obj.mcmc_type==1
            obj.sigma_mcmc = min(0.9999, obj.sigma_mcmc);
        end
        
        obj.Sigma = obj.Sigma_hat;
    elseif i > burnin
        obj.Sigma = obj.Sigma_hat;
    end
end

theta = fminsearch(@(x) -obj.logpi(x), theta);
theta = fminunc(@(x) -obj.logpi(x), theta);