function [Q, A, L] =  Q_RW1( dt, tau, y_index)
%%
%
%   RW2 precision operator
%    (x_t - x_{t-1}) \sim N(0, dt/tau )
%   x0 ~ N(0,\sigma1^2)
%   dt       - (double <1) the relative time step of dt_x= (dt * dt_y)
%   y_index  - (N x 1)     where y is observed
%
%
%%
N = length(y_index);
N_x = ceil(N/dt) + 1;
L = triu(sptoeplitz([-1 1 zeros(1,N_x  - 2 )]));
L =    L(1:end-1,:);
%L = (1 / dt^2) * sptoeplitz([-2 1 zeros(1,n )]); % second order operator
%L= sqrt(tau) *L(2:end-1,1:end);
Q =  tau/dt *( L' *  L ); % dt - due to the noise is N(0,1/dt)

A = speye(N_x);

Yobs = 2 + floor((0:(N-1))/dt);
A = A(Yobs,:);
A = A(y_index==1,:);