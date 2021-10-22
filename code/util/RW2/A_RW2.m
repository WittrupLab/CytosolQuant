function [Qt, At] =  A_RW2(Qx, dt, y_index)
%%
%   creating the observation matrix as well as
%   the subset of Qx needed to model y
%
%
%   Qx       - (n_x x n_x) preicision matrix
%   dt       - (double <1) the relative time step of dt_x= (dt * dt_y)
%   y_index  - (N x 1)     where y is observed
%
%
%
%%


N = length(y_index);
N_x = ceil(N/dt) + 2;
Qt = Qx(1:N_x, 1:N_x);

At = speye(N_x);

Yobs = 3 + floor((0:(N-1))/dt);
At = At(Yobs,:);
At = At(y_index==1,:);

