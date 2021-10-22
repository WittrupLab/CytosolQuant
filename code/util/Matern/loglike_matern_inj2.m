function lik = loglike_matern_inj2(theta, theta_inj, Y, inj, loc)
%%
%
%   Two matern random fields plus determinstic function
%   
%   theta - (7 x 1) parameter for noise and fields
%           (1)   mu      - mean 
%           (2-3) kappa   - log(range)
%           (4-5) sigma   - log(standard dev)
%           (6-7) nu      - log(differntiablility)
%           (8)   sigma_e - log(measurement error)
%   theta_inj (6 x 1) parameter for determistic function
%   D: 2019-01-19
%%



ind = ~isnan(Y);
[~, T] = size(Y);

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

lik = 0;


for i=1:T
    
  f_index    = find(ind(:,i));
  y_start  = min(f_index);
  y_end    = max(f_index);
  N = length(y_start:y_end);
  N_x = ceil(N/dt);
  Yobs =  1 + floor((0:(N-1))/dt);
  A = speye(N_x);
  A = A(Yobs,:);
  A = A(ind(y_start:y_end,i),:);  
  time     = (y_start + (0:(N_x-1))*dt)';
  [f, f_int] = injection_func(time, inj(i), t0, delta, beta, scale, dt);
  Ytb = Y(ind(:,i),i) - mu - A * f_int;
  Sigmai = Sigma(ind(:,i),ind(:,i));
  [R , p1] = chol(Sigmai);
  if(p1>0)
    lik = inf;
    return
  end
  v =  R'\Ytb;
  lik = lik - sum(log(diag(R))) - 0.5*(v'*v);
  
%     if(i==10)
%       
%        subplot(211)
%        plot(time, f)
%        subplot(212)
%        plot(time, f_int)
%        drawnow
%        
%     end
end
lik = (lik-0.5*size(R,1)*log(2*pi));