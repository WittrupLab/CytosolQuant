load('data/SIRNA_run_siGFP-1.mat')
sBIC = sum(exp(-BIC),2);
wBIC  = bsxfun(@rdivide, exp(-BIC), sBIC);
[m,idx] = max(wBIC,[],2);
index = sub2ind(size(R2),(1:size(R2,1))',idx);
plot( fitted_siRNA(index),R2_full(index),'.')
