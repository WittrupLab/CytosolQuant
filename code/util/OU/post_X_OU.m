function [X, E] = post_X_OU(theta, Y, dt)
%%
%   generating the posterior of X and noise
%   eps
%
%
%   theta - [mu, sigma, tau, sigma_eps]
%   mu        - (2 x 1) prior mean for x_0, x_1
%   sigma     - (2 x 1) prior cov for x_0, x_1
%   tau       - (2 x 1) precision scaling for RW2
%   sigma_eps - (1 x 1) measurement error
%   
%   Y     - (n   x n_obs) observations
%   dt       - (double <1) the relative time step of dt_x= (dt * dt_y)
%
%
%   X     - (n+2 x n_obs) latent processes
%   E     - (n+2 x n_obs) the noise driving the processes
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

X = NaN * zeros(size(Y));
E = cell(n_obs,4);

for i =1:n_obs
    index      = isnan(Y(:,i))==0;
    f_index    = find(index);
    y_start    = min(f_index);
    y_end      = max(f_index);
    Y_t        = Y(y_start:y_end,i);
   [Qx_t, At, mu_vec, D] =  Q_OU( dt, tau, lambda, mu, mu0, sigma0, index(y_start:y_end));
    
    
    y_t      = Y_t(isnan(Y_t)==0);
    Q_tilde  = (1/sigma_eps)^2*(At'*At) + Qx_t ;
    R_tilde  = chol(Q_tilde);
    v        = R_tilde'\(At' * (y_t/sigma_eps^2) + Qx_t * mu_vec);
    X_temp    = R_tilde \ v;
    
    X(index, i)  = At * X_temp;
    E{i,1} = D * X_temp;
    E{i,2} = y_start + (-2:(length(mu_vec)-3))*dt;
    E{i,3} = X_temp;
    E{i,4} = (X_temp(2:end)-X_temp(1:end-1))/dt;
end

