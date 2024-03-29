function bland_altman

load user_stats.mat

y = [user_stats.depth_abs_error];
x = [user_stats.comp_depth_real];
perc_90 = prctile(abs(y),90);
perc_50 = prctile(abs(y),50);
perc_10 = prctile(abs(y),10);

close all
figure(1)
subplot(2,1,1)
plot(x,y,'.')
hold on
plot([2 7],[0 0],'--r')
plot([2 7],[1 1]*perc_90,'Color',[1 1 1]*0.67,'LineStyle','-.')
plot([2 7],[1 1]*-perc_90,'Color',[1 1 1]*0.67,'LineStyle','-.')
plot([2 7],[1 1]*perc_50,'Color',[1 1 1]*0.5,'LineStyle','-.')
plot([2 7],[1 1]*-perc_50,'Color',[1 1 1]*0.5,'LineStyle','-.')
plot([2 7],[1 1]*perc_10,'Color',[1 1 1]*0.33,'LineStyle','-.')
plot([2 7],[1 1]*-perc_10,'Color',[1 1 1]*0.33,'LineStyle','-.')
title('Bland-Altman para el error de profundidad')
xlabel('Profundidad de compresión real (cm)')
ylabel('Error (cm)')
legend({'Error','Cero','Percentil 90','Percentil 50','Percentil 10'})

y = [user_stats.rate_abs_error];
x = [user_stats.comp_rate_real];
perc_90 = prctile(abs(y),90);
perc_50 = prctile(abs(y),50);
perc_10 = prctile(abs(y),10);

subplot(2,1,2)
plot(x,y,'.')
hold on
plot([60 140],[0 0],'--r')
plot([60 140],[1 1]*perc_90,'Color',[1 1 1]*0.67,'LineStyle','-.')
plot([60 140],[1 1]*-perc_90,'Color',[1 1 1]*0.67,'LineStyle','-.')
plot([60 140],[1 1]*perc_50,'Color',[1 1 1]*0.5,'LineStyle','-.')
plot([60 140],[1 1]*-perc_50,'Color',[1 1 1]*0.5,'LineStyle','-.')
plot([60 140],[1 1]*perc_10,'Color',[1 1 1]*0.33,'LineStyle','-.')
plot([60 140],[1 1]*-perc_10,'Color',[1 1 1]*0.33,'LineStyle','-.')
title('Bland-Altman para el error de frecuencia')
xlabel('Frecuencia de compresión real (cpm)')
ylabel('Error (cpm)')
legend({'Error','Cero','Percentil 90','Percentil 50','Percentil 10'})

end