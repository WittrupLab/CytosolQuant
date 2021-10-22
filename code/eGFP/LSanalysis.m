%%
%   Simple LS or abs to minize function
%   the decay function
%   D: 2019-11-15
%%

clearvars
close all
addpath(genpath('../util/'))

name = 'datasiGFP-1norm';%'datasiGFP-1norm'; % datasiGFP-2norm
load(name);
 theta_inj =[ 3
   -3
   -1
    3
    3
 ];
time = 500; % Event time
log_fit = 0;
dt = 1/10;
time_start = time -10;
time_end   = time + 400;
 %range of the siRNA values
 

 
%log distribution 
y = eGFP_traces;
if( log_fit)
    y(y<0) = 10^-4;
    y = log(y);
end
N = length(time_start:time_end);
N_x = ceil(N/dt);
t     = (time_start +  (0:(N_x-1))*dt)';
Yobs =  1 + floor((0:(N-1))/dt);
A = speye(N_x);
A = A(Yobs,:);
f = @(x)  optim_injection_func([time;x], t, A, y(time_start:time_end,:), siRNA_traces, log_fit);
theta = fminsearch(f, theta_inj);%

%theta = fminunc(f, theta);
theta = fminsearch(f, theta);%
theta = fminsearch(f, theta);%

save(['theta_',name,'.mat'], 'theta');
