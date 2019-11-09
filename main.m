clc;clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%LECTURA DATOS%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%GROUND TRUTH%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Nos devuelve una TABLA no matriz
M = cleanData('DATASET/anne1.csv');
%Separamos en arrays cada componente t, Ax, Ay, Az
%Columna 3 --> start Tiempo s?
%Columna 4 --> Depth mm?
t_real = M(:, 5).';
Depth_real = M(:, 3).';

%Pasamos a cm los mm
Depth_real = Depth_real * 10 ^ -1;

%%%%%%%%%%TEST DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Nos devuelve una TABLA no matriz
%Columna 1 --> Numero muestra
%Columna 2 --> Tiempo (ms)
%Columna 3,4,5 --> Acceleracion x, y, z (m/s^2)

N = cleanData('DATASET/precalib1_W.csv');

%Separamos en arrays cada componente t, Ax, Ay, Az
t = N(:, 1).';
Ax = N(:, 2).';
Ay = N(:, 3).';
Az = N(:, 4).';

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

%%Filtramos la señal de aceleracion para ELIMINAR la CONTINUA
A_filtrada = filtrado_senal(A, 1, 7, Fs, "high");
A_filtrada = filtrado_senal(A_filtrada, 40, 4, Fs, "low");

%Integral TRAPEZOIDAL
V = cumtrapz(t, A_filtrada);

%Filtramos la señal integrada (V) para ELIMINAR la CONTINUA
V_filtrada = filtrado_senal(V, 1, 7, Fs, "high");
V_filtrada = filtrado_senal(V_filtrada, 40, 4, Fs, "low");

%Pasamos a cm/s
V_filtrada = V_filtrada * 10 ^ 2;

%Volvemos a integrar para obtener la distancia
P = cumtrapz(t, V_filtrada);

%%Obtenemos los picos MAXIMOS y MINIMOS
[picos_max n_max] = findpeaks(P);
%Genio el de google
[picos_min n_min] = findpeaks(-P);

%Gestion de EXCEPCIONES
if length(n_max) > length(n_min)
    n_max = n_max(1:length(n_min));
    picos_max = picos_max(1:length(n_min));
end

if length(n_max) < length(n_min)
    n_min = n_min(1:length(n_max));
    picos_min = picos_min(1:length(n_max));
end

t_teorico = (n_max*1/Fs + n_min*1/Fs)/2;

%Y obtenemos finalmente la PROFUNDIDAD
Depth_teorica = picos_max + picos_min;

%Encuadrar Ventana (WINDOW) de trabajo
%%Limite INFERIOR
for i = 1:length(t_teorico)
    if t_teorico(1,i) > t_real(1,1)
        lim_inf = i;
        break;
    end
end

%%Limite SUPERIOR
W = 140;

%Calculamos el error de estimación
error = Depth_teorica(lim_inf:lim_inf + W - 1) - Depth_real(1:W);

%Calculamos el error cuadrático medio
err = mean((Depth_teorica(lim_inf:lim_inf + W - 1) - Depth_real(1:W)).^2);
disp(['El error cuadrático medio es: ',num2str(err)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%REPRESENTACION DATOS%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
subplot(3,1,1);
stem(t_teorico(lim_inf:lim_inf + W - 1), Depth_teorica(lim_inf:lim_inf + W - 1), '.');
title("Profundidad estimada")
xlabel("Instante tiempo (s)")
ylabel("Profundidad (cm)")
axis([0 max(t_teorico(lim_inf:lim_inf + W - 1)) + 2 0 8])

subplot(3,1,2);
stem(t_real(1:W), Depth_real(1:W), '.');
title("Profundidad real")
xlabel("Instante tiempo (s)")
ylabel("Profundidad (cm)")
axis([0 max(t_real(1:W)) + 2 0 8])

subplot(3,1,3);
stem(t_real(1:W), abs(error), '.');
title("Error")
xlabel("Compresión")
ylabel("Error (cm)")
axis([0 max(t_real(1:W)) + 2 0 8])

function dibujarP1(t, Ax, Ay, Az, A_filtrada, V, P);
    %Dibujamos cada ACELERACION por separado
    figure(1)

    subplot(6,1,1);
    plot(t,Ax);
    xlabel('Time (s)'); ylabel('Ax (m/s^2)');
    axis([t(1) t(length(t)) min(Ax) - 1 max(Ax) + 1]);

    subplot(6,1,2);
    plot(t,Ay);
    xlabel('Time (s)'); ylabel('Ay (m/s^2)');
    axis([t(1) t(length(t)) min(Ay) - 1 max(Ay) + 1])

    subplot(6,1,3);
    plot(t,Az);
    xlabel('Time (s)'); ylabel('Az (m/s^2)');
    axis([t(1) t(length(t)) min(Az) - 1 max(Az) + 1])

    subplot(6,1,4);
    plot(t,A_filtrada);
    xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
    axis([t(1) t(length(t)) min(A_filtrada) - 1 max(A_filtrada) + 1])

    subplot(6,1,5)
    plot(t,V);
    xlabel('Time (s)'); ylabel('Velocity (cm/s)');
    axis([t(1) t(length(t)) min(V) - 0.1 max(V) + 0.1])

    subplot(6,1,6)
    plot(t,P);
    xlabel('Time (s)'); ylabel('Distance (cm)');
    axis([t(1) t(length(t)) min(P) - 1 max(P) + 1])
end

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

function x = cleanData(nombre_fich)
    M = readtable(nombre_fich);
    N = table2array(M(:,2:width(M)));
    x = N;
end

