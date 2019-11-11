clc;clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%LECTURA DATOS%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%GROUND TRUTH%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Nos devuelve una TABLA no matriz
M = readtable('DATASET/anne2.csv');
%Separamos en arrays cada componente t, Ax, Ay, Az
%Columna 3 --> start Tiempo s?
%Columna 4 --> Depth mm?
t_real = ((M.startTime + M.endTime)/2)';
Depth_real = M.compDepth';

%Pasamos a cm los mm
Depth_real = Depth_real * 10 ^ -1;

%%%%%%%%%%TEST DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Nos devuelve una TABLA no matriz
%Columna 1 --> Numero muestra
%Columna 2 --> Tiempo (ms)
%Columna 3,4,5 --> Acceleracion x, y, z (m/s^2)

N = readtable('DATASET/precalib2_W.csv');

%Separamos en arrays cada componente t, Ax, Ay, Az
t = N.Time_Date';
Ax = N.Acc_X';
Ay = N.Acc_Y';
Az = N.Acc_Z';

%Pasamos a segundos los ms
t = t * 10 ^ -3;
%Hacemos que el eje empiece en 0
t = t - t(1);

%Obtenemos el periodo de muestra
Ts = t(2) - t(1);
%Obtenemos la frecuencia de muestreo
Fs = 1 / Ts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%PROCESADO DE DATOS%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Obtenemos la magnitud del vector
A = sqrt(Ax .^ 2 + Ay .^ 2 + Az .^ 2);

%%Filtramos la se�al de aceleracion para ELIMINAR la CONTINUA
A_filtrada = filtrado_senal(A, 1, 7, Fs, "high");
A_filtrada = filtrado_senal(A_filtrada, 40, 4, Fs, "low");

%Integral TRAPEZOIDAL
V = cumtrapz(t, A_filtrada);

%Filtramos la se�al integrada (V) para ELIMINAR la CONTINUA
V_filtrada = filtrado_senal(V, 1, 7, Fs, "high");
V_filtrada = filtrado_senal(V_filtrada, 40, 4, Fs, "low");

%Pasamos a cm/s
V_filtrada = V_filtrada * 10 ^ 2;

%Volvemos a integrar para obtener la distancia
P = cumtrapz(t, V_filtrada);

%%%QUITAMOS el TRANSITORIO --> Revisar
%Eliminamos un trozo de la se�al?
P = P(2*Fs:end);

%%Obtenemos los picos MAXIMOS y MINIMOS (Y) y la muestra en la que ocurre (X)
[picos_max n_max] = findpeaks(P);
%Genio el de google
[picos_min n_min] = findpeaks(-P);

%Gestion de EXCEPCIONES
%Si los arrays no son iguales, igualar dimensiones
if length(picos_max) > length(picos_min)
   picos_min = [picos_min, zeros(1,length(picos_max)-length(picos_min))];
   n_min = [n_min, zeros(1,length(n_max)-length(n_min))];
elseif length(picos_min) > length(picos_max)
   picos_max = [picos_max, zeros(1,length(picos_min)-length(picos_max))];
   n_max = [n_max, zeros(1,length(n_min)-length(n_max))];
end 

%Y obtenemos finalmente la PROFUNDIDAD
Depth_teorica = picos_max + picos_min;

t_teorico = (n_max*1/Fs + n_min*1/Fs)/2;

%Alinear se�al TEORICA y se�al REAL
dif = t_real(1,1) - t_teorico(1,1);
t_real = t_real - dif;

%Calculamos el error de estimaci�n
%error = Depth_teorica(lim_inf:lim_inf + W - 1) - Depth_real(1:W);

%Calculamos el error cuadr�tico medio
%err = mean((Depth_teorica(lim_inf:lim_inf + W - 1) - Depth_real(1:W)).^2);
%disp(['El error cuadr�tico medio es: ',num2str(err)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%REPRESENTACION GRAFICA%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
subplot(2,1,1);
stem(t_teorico, Depth_teorica, '.');
title("Profundidad estimada")
xlabel("Instante tiempo (s)")
ylabel("Profundidad (cm)")
axis([0 max(t_teorico) + 2 0 8])

subplot(2,1,2);
stem(t_real, Depth_real, '.');
title("Profundidad real")
xlabel("Instante tiempo (s)")
ylabel("Profundidad (cm)")
axis([0 max(t_real) + 2 0 8])

%subplot(3,1,3);
%stem(t_real, abs(error), '.');
%title("Error")
%xlabel("Compresi�n")
%ylabel("Error (cm)")
%axis([0 max(t_real) + 2 0 8])

function [P,f]= analisis_frecuencial(senal, Fs)
    %Respuesta en frecuencia
    A = senal.*hann(length(senal))';
    L = nextpow2(length(A));
    [X, f] = freqz(A, 1, 2^L, Fs);
    P = 1/2^L*abs(X).*2;
    f = f;
end

function senal_salida = filtrado_senal(senal_entrada, f_filtrado, orden, Fs, tipo)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%ELEGIR FILTRO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %tipo = "high" filtro paso-alto
    %tipo = "low" filtro paso-bajo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(tipo == "high")
        %%Definimos filtro IIR PASO-ALTO
        f0 = f_filtrado;
        [b a] = butter(orden, 2*f0/Fs, 'high');
    end
    
    if(tipo == "low")
        %%Definimos filtro PASO-BAJO
        f1 = f_filtrado;
        [b a] = butter(orden, 2*f1/Fs, 'low');
    end
    
    senal_salida = filtfilt(b,a,senal_entrada);
end

