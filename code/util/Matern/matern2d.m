function Sigma = matern2d(loc,kappa,sigma,nu)
%%
%	computing the matren covariance for a univariate random field
%
% INPUTS:
%  loc   - (n x 1) Locations to compute covariances at.
%  kappa - range parameter
%  nu    - smoothness parameter
%  sigma - std
%%


locs = ones(length(loc),2);
locs(:,1) = loc;
dist = squareform(pdist(locs));

Sigma = sigma(1)^2*materncorr(dist,kappa(1),nu(1));
