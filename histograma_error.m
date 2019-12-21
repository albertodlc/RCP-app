function histograma_error

load user_stats.mat

depth_abs_err = abs([user_stats.depth_abs_error]);
rate_abs_err = abs([user_stats.rate_abs_error]);

depth_rel_err = abs([user_stats.depth_rel_error]);
rate_rel_err = abs([user_stats.rate_rel_error]);

close all
figure(1)
subplot(2,1,1)
histogram(depth_abs_err,100);
title('Histograma del error absoluto en profundidad de compresi�n')
xlabel('|Error| (cm)')
ylabel('N�mero de casos')

subplot(2,1,2)
histogram(rate_abs_err,100);
title('Histograma del error absoluto en frecuencia de compresi�n')
xlabel('|Error| (cpm)')
ylabel('N�mero de casos')

figure(2)
subplot(2,1,1)
histogram(depth_rel_err,100);
title('Histograma del error relativo en profundidad de compresi�n')
xlabel('|Error| (%)')
ylabel('N�mero de casos')

subplot(2,1,2)
histogram(rate_rel_err,100);
title('Histograma del error absoluto en frecuencia de compresi�n')
xlabel('|Error| (%)')
ylabel('N�mero de casos')

end