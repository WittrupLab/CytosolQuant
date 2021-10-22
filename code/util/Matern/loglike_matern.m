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
kappa   = exp(theta(2));
sigma   = exp(theta(3));
nu      = min(exp(theta(4)),20);
sigma_e = exp(theta(5));
Sigma   = matern2d(loc',kappa,sigma,nu);
Sigma   = diag(length(Sigma))*sigma_e^2 + Sigma;


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