%%
%
% plotting dose response
%
% D: 2019-11-16
%%
close all
clearvars
addpath ../util/
addpath ../anders_code/boundedline/
addpath ../anders_code/Inpaint_nans/
time_start = 500;
time_vis = 600;
log_fit = 2;
load(sprintf('theta_vec_log%d.mat',log_fit));
n_sub = size(theta_vec,1);
dose_response = logspace(-3,3, 1000);



dt = 1/10;

N = length(time_start:time_vis);
N_x = ceil(N/dt);
t     = (time_start +  (0:(N_x-1))*dt)';
Relative_expression = zeros(length(dose_response), n_sub);
for i = 1:n_sub
    for ii=1:length(dose_response)
        [f,f_int] = injection_func(t, ...
                                   dose_response(ii),...
                                   time_start, ...
                                   exp(theta_vec(i,1)), ...
                                   exp(theta_vec(i,2:3)), ...
                                   exp(theta_vec(i,4:5)),...
                                   dt);
        Relative_expression(ii,i) = exp(f_int(end));
    end
end
Rel = quantile(Relative_expression,[0.025,0.5,0.975],2);
plot(log10(dose_response), Rel(:,2));
hold on
plot(log10(dose_response), Rel(:,1),'--');
plot(log10(dose_response), Rel(:,3),'--');
xlabel('dose log(siRNA)')
ylabel('Relative level')
title(sprintf('Time = %d',time_vis))
