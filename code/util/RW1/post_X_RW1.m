function [X, E] = post_X_RW1(theta, Y, dt)
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

mu          = theta(1);
sigma       = exp(theta(2));
tau         = exp(theta(3));
sigma_eps   = exp(theta(4));



[T_max, n_obs]  =size(Y);
mu = [mu; zeros(ceil(T_max/dt), 1)];
X = NaN * zeros(size(Y));
E = cell(n_obs,4);
for i =1:n_obs
    index    = isnan(Y(:,i))==0;
    f_index    = find(index);
    y_start  = min(f_index);
    y_end    = max(f_index);
    Y_t      = Y(y_start:y_end,i);
    
    [Qx_t, At, L] =  Q_RW1( dt, tau, index(y_start:y_end));
    n_t      = length(Qx_t);
    Q0 = sparse(n_t, n_t);
    Q0(1,1) = 1 /sigma(1)^2;
    y_t      = Y_t(isnan(Y_t)==0);
    mu_t     = mu(1:n_t);
    Q_tilde  = (1/sigma_eps)^2*(At'*At) + Qx_t + Q0;
    R_tilde  = chol(Q_tilde);
    v        = R_tilde' \ (At' * y_t/sigma_eps^2 + Q0*mu_t);
    X_temp    = R_tilde \ v;

    X(index, i)  = At * X_temp;
    E{i,1} = L * X_temp;
    E{i,2} = y_start + (-1:(n_t-3))*dt;
    E{i,3} = X_temp;
    E{i,4} = (X_temp(2:end)-X_temp(1:end-1))/dt;
end