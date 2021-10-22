function lik = loglike_inj_matern2(theta, theta_inj, Y, inj, loc)
%%
%   theta - (5 x 1)
%           mu      - mean
%           kappa   - range
%           sigma   - standard dev
%           nu      - differntiablility
%           sigma_e - measurement error
%
%%


T = length(inj);
ind = ~isnan(Y);

mu      = theta(1);
kappa   = exp(theta(2:3));
sigma   = exp(theta(4:5));
nu      = min(exp(theta(6:7)),20);
sigma_e = exp(theta(8));
Sigma   = matern2d(loc',kappa(1),sigma(1),nu(1));
Sigma2  = matern2d(loc',kappa(2),sigma(2),nu(2));
Sigma   = diag(length(Sigma))*sigma_e^2 + Sigma + Sigma2;


t0          = theta_inj(1);
delta       = exp(theta_inj(2));
beta        = exp(theta_inj(3:4));
scale       = exp(theta_inj(5:6));

dt=1/10;


f_index    = find(ind);
y_start  = min(f_index);
y_end    = max(f_index);
N = length(y_start:y_end);
N_x = ceil(N/dt);
Yobs =  1 + floor((0:(N-1))/dt);
A = speye(N_x);
A = A(Yobs,:);
A = A(ind(y_start:y_end),:);  
time     = (y_start + (0:(N_x-1))*dt)';
Sigmai = Sigma(ind, ind);
[R , ~] = chol(Sigmai);
lik = zeros(T,1);
for i=1:T
    
 
  [~, f_int] = injection_func(time, inj(i), t0, delta, beta, scale, dt);
  Ytb = Y(ind) - mu - A * f_int;
  
 
  v =  R'\Ytb;
  lik(i) =   - 0.5*(v'*v);
  
  
end