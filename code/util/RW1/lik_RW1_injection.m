function lik = lik_RW1_injection(theta, theta_inj, inj, Y, dt)
%%
%   loglikelhiood for RW1 model
%
%
%   theta - [mu, sigma, tau, sigma_eps]
%   mu        - (2 x 1) prior mean for x_0
%   sigma     - (2 x 1) prior cov for  x_0
%   tau       - (2 x 1) precision scaling for RW
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

mu          = theta(1);
sigma       = exp(theta(2));
tau         = exp(theta(3));
sigma_eps   = exp(theta(4));

t0          = theta_inj(1);
delta       = exp(theta_inj(2));
beta        = exp(theta_inj(3:4));
scale       = exp(theta_inj(5:6));

[T_max, n_obs]  =size(Y);
% T_max/dt gives the  number of discretization points, +2 is beacuse
% start two behind
mu = [mu; zeros(ceil(T_max/dt), 1)];
lik = 0;

for i =1:n_obs
    index    = isnan(Y(:,i))==0;
    f_index    = find(index);
    y_start  = min(f_index);
    y_end    = max(f_index);
    Y_t      = Y(y_start:y_end,i);
    [Qx_t, At] =  Q_RW1( dt, tau, index(y_start:y_end));
    
    n_t      = length(Qx_t);
    Q0 = sparse(n_t, n_t);
    Q0(1,1) = 1 /sigma(1)^2;
    
    time     = (y_start + (-1:(n_t-2))*dt)';
    [f, f_int] = injection_func(time, inj(i), t0, delta, beta, scale, dt);
     y_t      = Y_t(isnan(Y_t)==0) - At*f_int;
    mu_t     = mu(1:n_t);
    Q_tilde  = (1/sigma_eps)^2*(At'*At) + Qx_t + Q0;
    R_tilde  = chol(Q_tilde);
    v        = R_tilde'\(At' * (y_t/sigma_eps^2) + Q0*mu_t);
    lik      = lik + 0.5*(v'*v) - 0.5*(y_t'*y_t)/sigma_eps^2;
    lik      = lik - sum(log(diag(R_tilde)));
    
    R_x      = chol(Qx_t + Q0);
    lik      = lik + sum(log(diag(R_x)))  - length(y_t) * log(sigma_eps);
    lik      = lik - 0.5 * (mu_t' * Q0 * mu_t );
   if(i==10)
        figure(1)
       X_temp    = R_tilde \ v; 
       v_old = ( At' * ((y_t+At*f_int)/sigma_eps^2)+ Q0 * mu_t);
       X_old  = (R_tilde \ (R_tilde'\v_old)); 
       
       subplot(311)
       plot(time, X_temp) 
       hold on
       plot(time, X_old,'r') 
       plot(Y(:,i),'ko')
       hold off
       subplot(312)
       plot(time, f)
       subplot(313)
       plot(time, f_int)
       drawnow
       
    end
end
lik = - lik;