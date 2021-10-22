function [Q, A, mu_vec,  L0] =  Q_OU( dt, tau, lambda, mu, mu0, sigma0, y_index)
%%
%
%   dX_t = (\mu - X_t) dt + dW_t
%   OU precision operator
%    (x_t - x_{t-1})/dt + X_{t-1} dt\lambda \sim N( dt \lambda \mu ,1/(tau*dt) )
%   x0 ~ N(0,\sigma1^2)
%   x1 ~ N(0,\sigma2^2)
%   dt       - (double <1) the relative time step of dt_x= (dt * dt_y)
%   y_index  - (N x 1)     where y is observed
%
%
%%
N = length(y_index);
N_x = ceil(N/dt) + 1;
I = sparse(N_x, N_x);
I(2:end,1:end-1) = lambda * dt * speye(N_x-1);


D =   triu(sptoeplitz([-1 1 zeros(1,N_x  - 1 )]));
D=(D(1:end-1,2:end));
L = D + I; 
mu_vec_in  = mu * dt * lambda * ones(N_x, 1);
mu_vec_in(1) = mu0;
mu_vec  = L\mu_vec_in;


%L = (1 / dt^2) * sptoeplitz([-2 1 zeros(1,n )]); % second order operator
%L= sqrt(tau) *L(2:end-1,1:end);
Di = tau/dt * speye(N_x);
Di(1,1) = 1/sigma0^2;
Q = ( L' * Di *  L ); % dt - due to the noise is N(0,1/dt)

L0 = sqrt(tau/dt)* L(2:end,:);
%d_ = Q-L0'*L0*tau/dt;
%d_(1,1) = d_(1,1) - 1/sigma0^2;
%full(d_(1:4,1:4))

A = speye(N_x);

Yobs = 2 + floor((0:(N-1))/dt);
A = A(Yobs,:);
A = A(y_index==1,:);