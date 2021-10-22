function lik = loglike_matern2(theta, Y, loc)
%%
%   theta - (5 x 1)
%           mu      - mean
%           kappa   - range
%           sigma   - standard dev
%           nu      - differntiablility
%           sigma_e - measurement error
%
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


lik = 0;

for i=1:T
  Ytb = Y(ind(:,i),i) - mu;
  Sigmai = Sigma(ind(:,i),ind(:,i));
  [R , p1] = chol(Sigmai);
  if(p1>0)
    lik = inf;
    return
  end
  v =  R'\Ytb;
  lik = lik - sum(log(diag(R))) - 0.5*(v'*v);
end
lik = (lik-0.5*size(R,1)*log(2*pi));