%%
%   Simple LS or abs to minize function
%   the decay function
%   D: 2019-11-15
%%

clearvars
close all
addpath(genpath('../util/'))
n_sample = 1000;
load('datasiGFP-1norm');
 theta_inj =[ 3.7951
   -3.8475
   -1.4840
    3.2067
    3.8114
 ];
time = 500; % raden som motsvarar injection
log_fit = 0;
dt = 1/10;
time_start = time -10;
time_end   = time + 200;
 %range of the siRNA values
 

 
%log distribution 
y = eGFP_traces;
n_traj = size(y,2);
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
theta_vec = zeros(n_sample,length(theta_inj));
for j = 1:n_sample
    fprintf('iter  = %d\n',j);
    index = datasample((1:n_traj)',n_traj);
    f = @(x)  optim_injection_func([time;x], t, A, y(time_start:time_end,index), siRNA_traces(index), log_fit);
    theta = fminsearch(f, theta_inj);
    theta = fminsearch(f, theta);%
    theta_vec(j,:) = theta;
    disp(theta)
end

save(sprintf('theta_vec_log%d.mat',log_fit) , 'theta_vec');
