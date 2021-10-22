%%
%   pipeline for input frame to of siRNA input
%   using constant function, expontial decay, double expontial decay
%
% D: 2020-10-13

%
%%


% put the frames in frame_normalize_selected to zero mean
close all
clearvars 
saveRes = 0;
plot_flag = 0; %1 = plot only best
               %0 = plot all 

buttonPressFlag = 1; %0 = continuous run

repeat =  5; % how many time to repeat the runs
mcmc_sim   = 5000;
estimate_functions = [1, 2, 3];   % 1 - constant function , i.e. no expontial decay
                              % 2 - unit expontial decay
                              % 3 - double expontial decay




%%
% loading data and cleaning data
%%
addpath ../util/
addpath ../
load ('../../data/200928NonNorm/siGFP-1.mat')
siRNA_traces = collectedData.timeShift.maskedData.bcTraces.siRNA;
eGFP_traces = collectedData.timeShift.maskedData.bcTraces.eGFP;
normsiRNA_traces = NaN(size(siRNA_traces));

nTraces = size(siRNA_traces,2);

%Normalize to siRNA=0 during preceeding 10 frames
for i=1:nTraces
    normsiRNA_traces(1:end,i) = siRNA_traces(1:end,i)-(nanmean(siRNA_traces(490:499,i)));
end

%%
%setuping up values
%%  
R2 = zeros(nTraces, length(estimate_functions));
R2_full = zeros(nTraces, length(estimate_functions));
fitted_siRNA = nan(nTraces, length(estimate_functions));
BIC  = inf*ones(nTraces, length(estimate_functions));
start_pos=  nan(nTraces, length(estimate_functions));

% looping over the traces
for i=1:nTraces
    if all(not(isnan(normsiRNA_traces(500:515, i))))    %avoid an all nan column
        maxNonNaNFrame = lastNaNFree(normsiRNA_traces(501:end, i)); %to avoid nans
        [m, mi] = max(normsiRNA_traces(500:510, i)); %find max value in 10 frame window
        lastFrameToFit = min(mi+9, maxNonNaNFrame); %fit from max frame and 9 additional frame but avoid nans
        
        xdata =  (-10:15);
        xdata_plot = -10:0.1:15;
        xdata_anders = 0:15;
        ydata_anders = transpose(normsiRNA_traces((500+0):(500+15), i)); 
        ydata =  normsiRNA_traces(500+xdata, i);
        Sres = sum( (ydata_anders - nanmean(ydata_anders) ).^2,'omitnan') ;
        Sres_full = sum( (ydata - nanmean(ydata) ).^2,'omitnan') ;
        plot(xdata, ydata,'o');
        xlabel('t')
        ylabel('siRNA')
        hold on
        function_decay  = zeros(length(xdata_plot), 3);
        for r = 1:repeat
            for j = 1:length(estimate_functions)
                error_flag= true;
                while error_flag
                    
                    fprintf('*\n')
                    try
                        estimate_function = estimate_functions(j); 
                        n = sum(isnan(ydata_anders)==0);
                        if estimate_function==1 % Constant func
                            theta = noDecay_MCMCburnin(ydata, xdata, mcmc_sim);
                            resnorms = sum((ydata_anders(:) - noDecay(xdata_anders,theta)).^2,'omitnan');
                            sigma2   = mean((ydata_anders(:) - noDecay(xdata_anders,theta)).^2,'omitnan');
                            resnorms2 = sum((ydata(:) - noDecay(xdata,theta)).^2,'omitnan');
                            fitted_siRNA_temp = exp(theta(2));
                            function_decay_temp =  noDecay(xdata_plot,theta);
                            if plot_flag==0
                                plot(xdata_plot, function_decay_temp,'g')
                            end
                        elseif estimate_function==2 % univariate slope functions
                            theta = expDecay_MCMCburnin(ydata, xdata, mcmc_sim);
                            resnorms = sum((ydata_anders(:) - expDecay(xdata_anders,theta)).^2,'omitnan');
                            resnorms2 = sum((ydata(:) - expDecay(xdata,theta)).^2,'omitnan');
                            sigma2   = mean((ydata_anders(:) - expDecay(xdata_anders,theta)).^2,'omitnan');
                            fitted_siRNA_temp = exp( theta(2)) +  exp(theta(4));
                            function_decay_temp =   expDecay(xdata_plot,theta); 
                            if plot_flag==0
                                plot(xdata_plot, function_decay_temp,'r')
                            end
                        elseif estimate_function==3 % dubble expontial fit function
                            theta = expDecay2_MCMCburnin(ydata, xdata, mcmc_sim, 0);
                            resnorms = sum((ydata_anders(:) - expDecayDouble(xdata_anders,theta)).^2,'omitnan');
                            resnorms2 = sum((ydata(:) - expDecayDouble(xdata,theta)).^2,'omitnan');
                            sigma2   = mean((ydata_anders(:) - expDecayDouble(xdata_anders,theta)).^2,'omitnan');
                            fitted_siRNA_temp  = exp( theta(2)) +  exp(theta(4)) + exp( theta(8));
                            function_decay_temp = expDecayDouble(xdata_plot,theta);
                            if plot_flag==0
                                plot(xdata_plot, function_decay_temp,'b')
                            end
                        end
                        p = length(theta);
                        BIC_temp =  n*log(sigma2) + n + log(n)*p;
                        if(BIC_temp < BIC(i,j))
                            BIC(i,j) = BIC_temp;
                            R2(i,j)=1-resnorms/Sres ;
                            R2_full(i,j)=1-resnorms2/Sres_full ; 
                            start_pos(i,j) = exp(theta(1));
                            fitted_siRNA(i,j) = fitted_siRNA_temp;
                            
                            function_decay(:,j) = function_decay_temp;
                        end
                        error_flag=false;
                    catch e %e is an MException struct
                        fprintf(1,'The identifier was:\n%s',e.identifier);
                        fprintf(1,'There was an error! The message was:\n%s',e.message);
                        % more error handling...
                    end
                end
            end
            fitted_siRNA(i,:)
            %
        end
        
        selected = exp(-BIC(i,:))./sum(exp(-BIC(i,:)));
        
        title(sprintf(' R2= %.2f, BIC = %.2f, model probability = %.2f, %.2f, %.2f',max(R2(i,:)),min(BIC(i,:)), selected))
          
        if plot_flag==1
           plot(xdata_plot, function_decay(:,1) ,'g') 
           plot(xdata_plot, function_decay(:,2) ,'r') 
           plot(xdata_plot, function_decay(:,3) ,'b') 
       end
       % drawnow
       if buttonPressFlag==1
            waitforbuttonpress()
       end
       hold off
    else
        fprintf('all non\n')
    end
  
end
xlim([-10.5,15.5])
tightfig()
if(saveRes)
    save('data/SIRNA_run.mat', 'BIC', 'R2','R2_full', 'fitted_siRNA','start_pos');
end
