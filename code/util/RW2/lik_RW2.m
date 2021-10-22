function lik = lik_RW2(theta, Y, dt)
%%
%   loglikelhiood for RW2 model
%
%
%   theta - [mu, sigma, tau, sigma_eps]
%   mu        - (2 x 1) prior mean for x_0, x_1
%   sigma     - (2 x 1) prior cov for  x_0, x_1
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

mu          = theta(1:2);
sigma       = exp(theta(3:4));
tau         = exp(theta(5));
sigma_eps   = exp(theta(6));


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
    [Qx_t, At] =  Q_RW2( dt, tau, index(y_start:y_end));
    n_t      = length(Qx_t);
    Q0 = sparse(n_t, n_t);
    Q0(1,1) = 1 /sigma(1)^2;
    Q0(2,2) = 1 /sigma(2)^2;
    
    y_t      = Y_t(isnan(Y_t)==0);
    mu_t     = mu(1:n_t);
    Q_tilde  = (1/sigma_eps)^2*(At'*At) + Qx_t + Q0;
    R_tilde  = chol(Q_tilde);
    v        = R_tilde'\(At' * (y_t/sigma_eps^2) + Q0*mu_t);
    lik      = lik - sum(log(diag(R_tilde))) + (v'*v)/2;
    
    R_x      = chol(Qx_t + Q0);
    lik      = lik + sum(log(diag(R_x)))  - length(y_t) * log(sigma_eps);
    %lik      = lik -  log(sigma(1)) - log(sigma(2));
    lik      = lik - 0.5 * (mu_t' * Q0 * mu_t + (y_t'*y_t)/sigma_eps^2);
end
lik = - lik;