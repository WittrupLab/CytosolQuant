function lik = log_lik_inj_RW1(inj, Y, dt, theta, theta_inj)
%%
%   loglikelhiood for RW1 model
%   with injection effect on the graident
%
%
%   inj   - (m x 1) number of injection points
%
%   theta - [mu, log(sigma), log(tau), log(sigma_eps)]
%            mu        - (1 x 1) prior mean for x_0, 
%            sigma     - (1 x 1) prior cov for  x_0, 
%            tau       - (2 x 1) precision scaling for RW1
%            sigma_eps - (1 x 1) measurement error
%
%   theta_inj - [t_inj, log(delta), log(beta), log(scale) ]
%             t_inj      - (1 x 1) time for injection no effect prior
%             delta      - (1 x 1) t_inj+delta time where the effect is
%                                  largerst
%             beta       - (2 x 1) effect of the injection inj^beta
%             scale      - (1 x 1) dissparence of the injection  
%   
%   Y     - (n x 1) observations
%   dt    - (double < 1) discretization level of X, is actualy dt_x =
%                        dt*dt_y
%
%%



if nargin < 4
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

index    = isnan(Y)==0;
f_index    = find(index);
y_start  = min(f_index);
y_end    = max(f_index);
Y_t      = Y(y_start:y_end);
[Qx_t, At] =  Q_RW1( dt, tau, index(y_start:y_end));
n_t      = length(Qx_t);
Q0 = sparse(n_t, n_t);
Q0(1,1) = 1 /sigma(1)^2;

mu = [mu; zeros(ceil(length(Y)/dt), 1)];
mu_t     = mu(1:n_t);
Q_tilde  = (1/sigma_eps)^2*(At'*At) + Qx_t + Q0;
R_tilde  = chol(Q_tilde);

time     = (y_start + (-1:(n_t-2))*dt)';
lik = zeros(length(inj),1);
%figure(20)
figure(1)

for i =1:length(inj)
    [f, f_int] = injection_func(time, inj(i), t0, delta, beta, scale, dt);
    y_t      = Y_t(isnan(Y_t)==0) - At*f_int;
    v        = R_tilde'\(At' * (y_t/sigma_eps^2) + Q0 * mu_t);
    lik(i)      =  -  0.5*(y_t'*y_t)/sigma_eps^2 + 0.5*(v'*v);


%     subplot(311)
%     plot(time, R_tilde\v);
%     hold on
%     subplot(312)
%     plot(time, f) 
%     hold on
%     subplot(313)
%     plot(time,f_int)
%     hold on
end