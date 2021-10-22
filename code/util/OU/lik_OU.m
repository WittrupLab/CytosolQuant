function lik = lik_OU(theta, Y, dt)
%%
%   loglikelhiood for RW2 model
%
%
%   theta - [mu, sigma, tau, sigma_eps]
%   mu        - (2 x 1) prior mean for x_0, x_1
%   sigma     - (2 x 1) prior cov for x_0, x_1
%   tau       - (2 x 1) precision scaling for RW2
%   sigma_eps - (1 x 1) measurement error
%   
%   Y     - (n x n_obs) observations
%   dt    - (double < 1) discretization level of X, is actualy dt_x =
%                        dt*dt_y
%
%%

if nargin < 3
    dt = 1;
end

mu0          = theta(1);
mu           = theta(2);
sigma0      = exp(theta(3));
lambda      = exp(theta(4));
tau         = exp(theta(5));
sigma_eps   = exp(theta(6));


[T_max, n_obs]  =size(Y);
% T_max/dt gives the  number of discretization points, +2 is beacuse
% start two behind
lik = 0;
lik_v = zeros(n_obs,1);
for i =1:n_obs
    index      = isnan(Y(:,i))==0;
    f_index    = find(index);
    y_start    = min(f_index);
    y_end      = max(f_index);
    Y_t        = Y(y_start:y_end,i);
    [Qx_t, At, mu_vec, L0] =  Q_OU( dt, tau, lambda, mu, mu0, sigma0, index(y_start:y_end));
    y_t      = Y_t(isnan(Y_t)==0) - At*mu_vec;
    Q_tilde  = (1/sigma_eps)^2*(At'*At) + Qx_t;
    R_tilde  = chol(Q_tilde);
    v        = R_tilde'\(At' * (y_t/sigma_eps^2) );
    lik_v(i) = lik;
    lik      = lik - sum(log(diag(R_tilde))) + (v'*v)/2;
    
    [R_x,p]      = chol(Qx_t);
    if p==1
        lik = inf;
       return; 
    end
    %n_x = length(Qx_t);
    lik      = lik + sum(log(diag(R_x))) - length(y_t) * log(sigma_eps);
    %lik      = lik - log(sigma0) + (n_x-1)*log(tau)/2 - length(y_t) * log(sigma_eps);
    lik      = lik - 0.5 *  (y_t'*y_t)/sigma_eps^2;
    lik_v(i) = lik -lik_v(i);
         figure(i)
         subplot(211)
        step= y_start + (-1:(length(mu_vec)-2))*dt;
        plot(step, (R_tilde\v))
        hold on
        Y_T = Y(:,i);
        Y_T(isnan(Y_T)==0) =Y_T(isnan(Y_T)==0) - At * mu_vec;
        plot(Y_T,'ro')
        hold off
        %plot(step, mu_vec,'g')
        subplot(212)
         plot(step,sqrt(diag(full(inv(Qx_t)))))
        drawnow
         
%        
end
lik = - lik;