
clearvars
addpath boundedline/
addpath Inpaint_nans/
load('../siRNA/data/SIRNA_run_siGFP-1.mat'); %load release quantifications
run('plotAllMeansNonNorm.m'); %loads inactive siRNA correction

saveMovie = 0; %set to 1 to save images of scatter movie

load ('../../data/200928NonNorm/siGFP-1.mat') %load siRNA and eGFP traces

% set indices based on lowest BIC from siRNA release modelling
sBIC = sum(exp(-BIC),2);
wBIC  = bsxfun(@rdivide, exp(-BIC), sBIC);
[m,idx] = max(wBIC,[],2);
index = sub2ind(size(R2),(1:size(R2,1))',idx);



%extract masked data from tracked cells
siRNA_traces = collectedData.timeShift.maskedData.bcTraces.siRNA;
eGFP_traces = collectedData.timeShift.maskedData.bcTraces.eGFP;

timeshiftedsiRNA = siRNA_traces;
trimmedeGFP = eGFP_traces;


releaseFluo = fitted_siRNA(index);
nEvents = length(releaseFluo);

%% Remove cells with poor siRNA quant model fit
poorQ = R2_full(index)<0.75;
trimmedeGFP(1:end, poorQ)=NaN;
nonQ = isnan(releaseFluo);
trimmedeGFP(1:end, nonQ)=NaN;


%% Remove cells with missed events 


minVal = min(trimmedeGFP(500:510,1:end));
missedEventCells=minVal < 0.001;
trimmedeGFP(1:end, missedEventCells)=NaN;



[sortedReleaseFluo, I] = sort(releaseFluo);
releaseSortedeGFP = trimmedeGFP(1:end,I);

[m, eventOffset] = max(timeshiftedsiRNA(500:510, :), [], 1);
eventOffset= eventOffset-1;
for i=1:nEvents   
       timeshiftedsiRNA(1:(end-eventOffset(i)),i) = timeshiftedsiRNA((1+eventOffset(i)):end,i);
end

releaseSortedsiRNA = timeshiftedsiRNA(1:end,I);

sorted_R2 = R2_full(index);

sorted_R2 = sorted_R2(I);




%% Correct inertial eGFP increase (cell growth)
inertiaCorrectedeGFP = releaseSortedeGFP./siLucCorr;


%% Trim away all NaN Traces

remainingTraces = not(isnan(inertiaCorrectedeGFP(500, 1:end)));
eGFPAnalysis = inertiaCorrectedeGFP(1:end, remainingTraces);
releaseFluoAnalysis = sortedReleaseFluo(remainingTraces);
siRNAAnalysis = releaseSortedsiRNA(1:end, remainingTraces);
eGFPRawAnalysis = eGFP_traces(1:end, I);
eGFPRawAnalysis = eGFPRawAnalysis(1:end, remainingTraces);
R2_analysis = sorted_R2(remainingTraces);

%%
nPools = 4;
nEvents=sum(remainingTraces);
poolsize=nEvents/nPools;


figure;
hold on
colororder = get(gca, 'colororder');
xlabel('Time (min)')
ylabel('Relative d1eGFP expression')
xlim([-50, 600])
ylim([0, 1.2])


releasePools = NaN(round(poolsize+0.5), nPools);


for i=1:nPools
    tempPool = eGFPAnalysis(1:end, ((i-1)*round(poolsize)+1):round(i*poolsize));
    tempPool = tempPool./nanmean(tempPool(500,1:end));
    means = nanmean(tempPool, 2);
    stdev = nanstd(transpose(tempPool));
    nTraces = transpose(sum(not(isnan(tempPool)), 2));
    CI80 = 1.28*stdev./sqrt(nTraces);
    boundedline(-499*5:5:500*5, means, CI80, 'alpha', 'cmap', colororder(mod(i-1,7)+1, 1:end), 'transparency', 0.1);
    median(releaseFluoAnalysis(((i-1)*round(poolsize)+1):round(i*poolsize)))
    releasePools(1:size(tempPool,2),i) = releaseFluoAnalysis(((i-1)*round(poolsize)+1):round(i*poolsize));
end 

hold off


%%
figure;
hold on

for i=1:nPools
    
    means = nanmedian(siRNAAnalysis(1:end, ((i-1)*round(poolsize)+1):round(i*poolsize)), 2);
    stdev = nanstd(transpose(siRNAAnalysis(1:end, ((i-1)*round(poolsize)+1):round(i*poolsize))));
    nTraces = sum(not(isnan(transpose(siRNAAnalysis(1:end, ((i-1)*round(poolsize)+1):round(i*poolsize))))), 1);
    CI80 = 1.28*stdev./sqrt(nTraces);
    boundedline(-499*5:5:500*5, means, CI80, 'alpha', 'cmap', colororder(i, 1:end), 'transparency', 0.1);
    
end
xlabel('Time (min)')
ylabel('[siRNA] (nM)')
xlim([-50, 150])
ylim([-1, 40])
hold off

% 
normeGFPAnalysis = eGFPAnalysis./eGFPAnalysis(500, 1:end);
for i=0:100
    figure(15);scatter(log10(releaseFluoAnalysis), normeGFPAnalysis(500+i, 1:end), R2_analysis.*50)
    xlabel('[siRNA] (log(nM))')
    ylabel('d1eGFP expression')
    ylim([0, 2])
    xlim([-1, 2.5])
    filename= sprintf('./Figures/ScatterMovie/a3343scatter%d.png', i);
    if saveMovie
        saveas(gcf, filename)
    end
end
