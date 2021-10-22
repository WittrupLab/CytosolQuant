%%%
% visualsing the functions
%
%%%


close all
clearvars
addpath ../QuantilePlots/boundedline/
addpath ../QuantilePlots/Inpaint_nans/

addpath(genpath('../util/'))
name = 'datasiGFP-1norm';  %datasiGFP-1norm datasiGFP-2norm
load(name);
load(['theta_',name,'.mat']);
%%
logscale = 0;
nPools = 5;
nEvents= size(eGFP_traces,2);
poolsize=nEvents/nPools;
average = 1;


theta_inj =[
    500  
    theta
 ];

time = theta_inj(1);
time_start = time -10;
time_end   = time + 120;

figure(nPools);

hold on

sim = 1001;
dt = 1/10;

N = length(time_start:time_end);
N_x = ceil(N/dt);
t     = (time_start +  (0:(N_x-1))*dt)';
colororder = get(gca, 'colororder');
if(logscale)
    eGFP_traces = log(eGFP_traces);
end
for i=1:nPools
    index = ((i-1)*round(poolsize)+1):round(i*poolsize);
    sI    =siRNA_traces(index);
    tempPool = eGFP_traces(time_start:time_end, index);
    stats = bootstrp(sim,@(x) nanmedian(x,1),tempPool')';
    means =nanmean(stats,2);
    EMP80 = quantile(stats,[0.1,0.9],2);
    rem = sum(isnan(tempPool)==0,2) < 5;
    EMP80(rem,:) = nan;
    means(rem) = nan;
    EMP80 = abs(bsxfun(@minus, EMP80 , means));
    %boundedline(time_start:time_end, means, EMP80, 'alpha', 'cmap', colororder(mod(i-1,7)+1, 1:end));
    median(siRNA_traces(((i-1)*round(poolsize)+1):round(i*poolsize)))
    if(average)
       [f,f_int] = injection_func(t, ...
                               median(sI),...
                               theta_inj(1), ...
                               exp(theta_inj(2)), ...
                               exp(theta_inj(3:4)), ...
                               exp(theta_inj(5:6)),...
                               dt);
        F = f_int;
    else
        F = zeros(length(t),length(sI)); 
        for j=1:length(sI)
            [f,f_int] = injection_func(t, ...
                                   sI(j),...
                                   theta_inj(1), ...
                                   exp(theta_inj(2)), ...
                                   exp(theta_inj(3:4)), ...
                                   exp(theta_inj(5:6)),...
                                   dt);
            F(:,j) = f_int;                  
        end
    end
   A = speye(length(t));
   if(logscale)
    plot(A*t,median(A*F,2),'color', colororder(mod(i-1,7)+1, 1:end));    
   else
    plot(((A*t)-500).*5,median(exp(A*F),2),'color', colororder(mod(i-1,7)+1, 1:end));    
   end
end 

hold off

for i=1:nPools
    index = ((i-1)*round(poolsize)+1):round(i*poolsize);
    tempPool = eGFP_traces(time_start:time_end, index);
     sI    =siRNA_traces(index);
figure(nPools+i)

[f,f_int] = injection_func(t, ...
                           median(siRNA_traces(index)),...
                           theta_inj(1), ...
                           exp(theta_inj(2)), ...
                           exp(theta_inj(3:4)), ...
                           exp(theta_inj(5:6)),...
                           dt);
   A = speye(length(t));
   plot(time_start:time_end,tempPool);
   hold on
    if(average)
       [f,f_int] = injection_func(t, ...
                               median(sI),...
                               theta_inj(1), ...
                               exp(theta_inj(2)), ...
                               exp(theta_inj(3:4)), ...
                               exp(theta_inj(5:6)),...
                               dt);
        F = f_int;
    else
        F = zeros(length(t),length(sI)); 
        for j=1:length(sI)
            [f,f_int] = injection_func(t, ...
                                   sI(j),...
                                   theta_inj(1), ...
                                   exp(theta_inj(2)), ...
                                   exp(theta_inj(3:4)), ...
                                   exp(theta_inj(5:6)),...
                                   dt);
            F(:,j) = f_int;                  
        end
    end

    if(logscale)
     plot((A*t),median((A*F),2),'color', colororder(mod(i-1,7)+1, 1:end),'lineWidth',4);
   else
    plot(A*t,median(exp(A*F),2),'color', colororder(mod(i-1,7)+1, 1:end),'lineWidth',4); 
   end
  
end
