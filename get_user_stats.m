clc;clear;close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%LECTURA DATOS%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%GROUND TRUTH%%%%%%%%%%%%%%%%%%%%%%%%%%%%

user_stats = struct();

for num = 1:15

%%Especificar desfase
% desf = input('Introduce desfase: '); %%Introducir desfase manualmente
load desf % Carga array con el desfase calculado para cada dataset

%%%Nos devuelve una TABLA no matriz
M = readtable(strcat('DATASET/anne',num2str(num),'.csv'));
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

N = readtable(strcat('DATASET/precalib',num2str(num),'_W.csv'));

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

%%Añadimos 2s de señal antes y después, para evitar el transitorio
n_trans = round(2*Fs);
A = [ones(1,n_trans)*A(1) A ones(1,n_trans)*A(end)];

%%Filtramos la señal de aceleracion para ELIMINAR la CONTINUA
A_filtrada = filtrado_senal(A, 1, 7, Fs, "high");
A_filtrada = filtrado_senal(A_filtrada, 40, 4, Fs, "low");

%%Eliminamos la señal añadida
A_filtrada = A_filtrada(n_trans+1:end-n_trans);

%Integral TRAPEZOIDAL
V = cumtrapz(t, A_filtrada);

%%Ponemos señal adicional para evitar transitorios
V = [ones(1,n_trans)*V(1) V ones(1,n_trans)*V(end)];

%Filtramos la señal integrada (V) para ELIMINAR la CONTINUA
V_filtrada = filtrado_senal(V, 1, 7, Fs, "high");
V_filtrada = filtrado_senal(V_filtrada, 40, 4, Fs, "low");

%Eliminamos señal adicional
V_filtrada = V_filtrada(n_trans+1:end-n_trans);

%Pasamos a cm/s
V_filtrada = V_filtrada * 10 ^ 2;

%Volvemos a integrar para obtener la distancia
P = cumtrapz(t, V_filtrada);

%%%QUITAMOS el TRANSITORIO --> Revisar
%Eliminamos un trozo de la señal?
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
if t_teorico(end) < t_teorico(end-1) %Cuando los picos no son pares, hay error en el último valor
   t_teorico = t_teorico(1:end-1); %Si hay error, se elimina
   Depth_teorica = Depth_teorica(1:end-1);
end

%Corregir desfase entre señal teórica y real
if desf(num) < 0
    t_teorico = t_teorico - desf(num);
    
    Depth_real = Depth_real(t_real>=t_teorico(1)); %Se elimina la parte de la señal real anterior al comienzo de la teórica
    t_real = t_real(t_real>=t_teorico(1));
    t_real = t_real - t_teorico(1); %La señal teórica comienza en t = 0
    t_teorico = t_teorico - t_teorico(1);
elseif desf(num) > 0
    t_real = t_real + desf(num);
    
    Depth_teorica = Depth_teorica(t_teorico>=t_real(1)); %Se elimina la parte de la señal teórica anterior al comienzo de la real
    t_teorico = t_teorico(t_teorico>=t_real(1));
    t_teorico = t_teorico - t_real(1); %La señal real comienza en t = 0
    t_real = t_real - t_real(1);
end

%%Igualar duración en segundos de la señal teórica y la real
if t_real(end) > t_teorico(end)
    Depth_real = Depth_real(t_real<=t_teorico(end));
    t_real = t_real(t_real<=t_teorico(end));
elseif t_teorico(end) > t_real(end)
    Depth_teorica = Depth_teorica(t_teorico<=t_real(end));
    t_teorico = t_teorico(t_teorico<=t_real(end));
end

t_window = 5;
step = 2.5;
[comp_rate_teorica, comp_depth_teorica, comp_t_teorica] = analisis_senal(Depth_teorica,t_teorico,t_window,step);
[comp_rate_real, comp_depth_real, comp_t_real] = analisis_senal(Depth_real,t_real,t_window,step);

if length(comp_t_real) > length(comp_t_teorica)
    comp_t_real = comp_t_real(1:end-1);
    comp_depth_real = comp_depth_real(1:end-1);
    comp_rate_real = comp_rate_real(1:end-1);
elseif length(comp_t_real) < length(comp_t_teorica)
    comp_t_teorica = comp_t_teorica(1:end-1);
    comp_depth_teorica = comp_depth_teorica(1:end-1);
    comp_rate_teorica = comp_rate_teorica(1:end-1);
end

rate_abs_error = comp_rate_teorica - comp_rate_real;
rate_rel_error = (comp_rate_teorica - comp_rate_real)*100./comp_rate_real;

depth_abs_error = comp_depth_teorica - comp_depth_real;
depth_rel_error = (comp_depth_teorica - comp_depth_real)*100./comp_depth_real;

user_stats(num).comp_rate_teorica = comp_rate_teorica;
user_stats(num).comp_rate_real = comp_rate_real;
user_stats(num).comp_depth_teorica = comp_depth_teorica;
user_stats(num).comp_depth_real = comp_depth_real;
user_stats(num).comp_t = comp_t_real;
user_stats(num).rate_abs_error = rate_abs_error;
user_stats(num).rate_rel_error = rate_rel_error;
user_stats(num).depth_abs_error = depth_abs_error;
user_stats(num).depth_rel_error = depth_rel_error;
end

save user_stats user_stats

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

function [comp_rate, comp_depth,comp_t] = analisis_senal(depth,t,t_window,step)
% Extraer valores medios de frecuencia y profundidad de compresión
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% depth: Array de profundidad de compresiones
% t: Tiempo de array de compresiones
% t_window: Tiempo de ventana
% step: Tiempo transcurrido entre cada feedback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Teniendo en cuenta que t_window > step
% El primer feedback se da cuando se tiene una señal de duración = t_window
% Se termina de analizar el tiempo que queda hasta depth(end) < step

    pasos = floor((t(end)-t_window)/step) + 1;

    comp_rate = zeros(1,pasos);
    comp_depth = zeros(1,pasos);
    comp_t = linspace(t_window,(pasos-1)*step+t_window,pasos);
      
    for i = 0:pasos-1
        depth_aux = depth(t>=i*step & t<(i*step + t_window));
        t_aux = t(t>=i*step & t<(i*step + t_window));
        
        comp_rate(i+1) = 60./median(diff(t_aux)); % En compresiones por minuto
        comp_depth(i+1) = median(depth_aux);
    end
end