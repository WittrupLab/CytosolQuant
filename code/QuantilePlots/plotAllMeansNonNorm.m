%Script to plot the average eGFP expression of all event cells for siGFP-1
%(potent siRNA), siGFP-2 (less potent siRNA) and siLuc (inactive siRNA). A
%moving average of siLuc treated cells (siLucCorr) is also calculated to compensate for
%baseline drift (initial apparent bump in eGFP signal) and corrected data
%is plotted.



load('../../data/200928NonNorm/siGFP-2.mat');
siGFP2=collectedData.timeShift.maskedData.bcTraces.eGFP;
load('../../data/200928NonNorm/siGFP-1.mat');
siGFP1=collectedData.timeShift.maskedData.bcTraces.eGFP;
load('../../data/200928NonNorm/siLuc.mat');
siLuc=collectedData.timeShift.maskedData.bcTraces.eGFP;
siLucCorr = movmean(nanmean(siLuc, 2), [5,5]);
siLucCorr = siLucCorr./siLucCorr(500);
siLucCorr(1:470)=siLucCorr(470);
siLucCorr(560:end)=siLucCorr(560);

siGFP1c=siGFP1./siLucCorr;
siGFP2c=siGFP2./siLucCorr;
siLucc=siLuc./siLucCorr;



figure
hold on
plot(nanmean(siGFP1, 2)./nanmean(siGFP1(500,:)))
plot(nanmean(siGFP2, 2)./nanmean(siGFP2(500,:)))
plot(nanmean(siLuc, 2)./nanmean(siLuc(500,:)))
plot(siLucCorr)
axis([470, 600, 0, 1.5])

figure
hold on
plot(nanmean(siGFP1, 2)./nanmean(siGFP1(500,:))./siLucCorr)
plot(nanmean(siGFP2, 2)./nanmean(siGFP2(500,:))./siLucCorr)
plot(nanmean(siLuc, 2)./nanmean(siLuc(500,:))./siLucCorr)
axis([470, 600, 0, 1.5])

figure(10);
hold on
colororder = get(gca, 'colororder');
colororder=colororder.*0.5; %darken
xlabel('Time (min)')
ylabel('Relative d1eGFP expression')
xlim([-80, 500])
ylim([0, 1.4])

allTraces={siGFP1, siGFP2, siLuc};

for i=1:3
    traces=allTraces{i};
    tracesNorm=traces./nanmean(traces(500,:));
    tracesStd=nanstd(tracesNorm, 0, 2);
    tracesObs=sum(not(isnan(tracesNorm)),2);
    CI95=(1.96*tracesStd)./sqrt(tracesObs);
    CI80=(1.28*tracesStd)./sqrt(tracesObs);
    boundedline(-499*5:5:500*5, nanmean(tracesNorm, 2), CI80, 'alpha', 'cmap', colororder(mod(i-1,7)+1, 1:end));
end

plot(-499*5:5:500*5,siLucCorr, 'k')
%legend('80%CI','siGFP-1', '80%CI', 'siGFP-2', '80%CI', 'siLuc', 'Correction')

figure(11);
hold on
colororder = get(gca, 'colororder');
colororder=colororder.*0.5; %darken
xlabel('Time (min)')
ylabel('Relative d1eGFP expression')
xlim([-80, 500])
ylim([0, 1.4])

allTraces={siGFP1c, siGFP2c, siLucc};

for i=1:3
    traces=allTraces{i};
    tracesNorm=traces./nanmean(traces(500,:));
    tracesStd=nanstd(tracesNorm, 0, 2);
    tracesObs=sum(not(isnan(tracesNorm)),2);
    CI95=(1.96*tracesStd)./sqrt(tracesObs);
    CI80=(1.28*tracesStd)./sqrt(tracesObs);
    boundedline(-499*5:5:500*5, nanmean(tracesNorm, 2), CI80, 'alpha', 'cmap', colororder(mod(i-1,7)+1, 1:end));
end
%legend('80%CI','siGFP-1', '80%CI', 'siGFP-2', '80%CI', 'siLuc')