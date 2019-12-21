function crear_boxplot

load user_stats.mat

close all
figure(1)

depth_rel_error = [user_stats.depth_rel_error];
rate_rel_error = [user_stats.rate_rel_error];
dataset_length = zeros(1,15);
g = cell(1,15);

%% A la función boxplot se le pasa un array con datos y otro array con etiquetas (strings)
%% con el nombre de la clase (en este caso dataset) el que pertenece cada dato
for i = 1:15
    dataset_length(i) = length(user_stats(i).depth_rel_error);
    g{i} = repmat({strcat('Dataset ',num2str(i))},1,dataset_length(i)); % Crea un array de celdas con el string 'Dataset #' en cada una
end

g_t = [g{1},g{2},g{3},g{4},g{5},g{6},g{7},g{8},g{9},g{10},g{11},g{12},g{13},g{14},g{15}]; % Todas las celdas con etiquetas puestas en fila

close all
figure(1)
boxplot(depth_rel_error,g_t)
title('Error relativo de profundidad')
ylabel('Error relativo (%)')

figure(2)
boxplot(rate_rel_error,g_t)
title('Error relativo de frecuencia')
ylabel('Error relativo (%)')
end

