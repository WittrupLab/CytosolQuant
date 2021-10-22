function [sigma, acc_sum]  = AMCMC_RR(sigma, acc_sum, iter, batch, accepantce_rate, delta_min, delta_rate)
%%
% adapting sigma from Roberts a Rosenthal
%
%   sigma           - (double>0)  acceptance prob
%   acc_sum         - (acc_sum)   number of accapted samples in batch samples
%   iter            - (int)       iteration number in the MCMC algorithm
%   batch           - (int d:50)  number samples on which to test acceptance rate
%   accepantce_rate - (double)    target acceptance
%   delta_min       - (double)    exp(delta_min) is largest multplicative
%                                 update
%   delta_rate      - (double)   decrease of adapation is
%                                   iter^(-delta_rate)
%
%
%  update: 2018-07-21
%%

if nargin < 4
    batch = 50;
end
if nargin < 5
    accepantce_rate = 0.234;
end
if nargin < 6
    delta_min =0.1;
end
if nargin < 7
    delta_rate = 1/3;
end

if mod(iter , batch) == 0
    delta = min(delta_min, (iter/batch)^(-delta_rate));
    if acc_sum/batch > accepantce_rate
       sigma = sigma * exp( delta);
    else
        sigma = sigma / exp( delta);
    end
    
    acc_sum = 0;  
end