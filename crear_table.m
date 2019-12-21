function crear_table

load user_stats.mat

Dataset = [1:15]';
Mediana_error_absoluto_profundidad = zeros(15,1);
Mediana_error_relativo_profundidad = zeros(15,1);
Mediana_error_absoluto_frecuencia = zeros(15,1);
Mediana_error_relativo_frecuencia = zeros(15,1);

for i = 1:15
    Mediana_error_absoluto_profundidad(i,1) = median(abs(user_stats(i).depth_abs_error));
    Mediana_error_relativo_profundidad(i,1) = median(abs(user_stats(i).depth_rel_error));
    Mediana_error_absoluto_frecuencia(i,1) = nanmedian(abs(user_stats(i).rate_abs_error)); % Para evitar un valor NaN que hay en Dataset3
    Mediana_error_relativo_frecuencia(i,1) = nanmedian(abs(user_stats(i).rate_rel_error));
end

table(Dataset,Mediana_error_absoluto_profundidad,Mediana_error_relativo_profundidad,...
    Mediana_error_absoluto_frecuencia,Mediana_error_relativo_frecuencia)

end

