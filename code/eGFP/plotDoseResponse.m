%%
%
% plotting dose response
%
% D: 2019-11-16
%%
close all
clearvars
time_vis = 600;
theta_inj =[
    500     
    3.6764
   -4.4466
   -0.9622
    2.9835
   25.1278
 ];
theta_inj =[
    500   
    2.9816
   -3.9106
   -1.4824
   -2.9131
    4.5601
 ];
dose_response = logspace(-6,4, 10000);



dt = 1/10;

N = length(theta_inj(1):time_vis);
N_x = ceil(N/dt);
t     = (theta_inj(1) +  (0:(N_x-1))*dt)';
Relative_expression = zeros(length(dose_response), 1);
for i=1:length(dose_response)
    [f,f_int] = injection_func(t, ...
                               dose_response(i),...
                               theta_inj(1), ...
                               exp(theta_inj(2)), ...
                               exp(theta_inj(3:4)), ...
                               exp(theta_inj(5:6)),...
                               dt);
    Relative_expression(i) = exp(f_int(end));  
end
plot(dose_response, Relative_expression);
xlabel('dose siRNA')
ylabel('Relative level')
title(sprintf('Time = %d',time_vis))
