clc;clear;
%%%Nos devuelve una TABLA no matriz
%Columna 1 --> Numero muestra
%Columna 2 --> Tiempo ms?
%Columna 3,4,5 --> Acceleracion x, y, z
M = readtable('precalib1_W.csv');
%Convertimos a array
N = table2array(M(:,2:5));
%Separamos en arrays cada componente t, Ax, Ay, Az
t = N(:,1).';
Ax = N(:,2).';
Ay = N(:,3).';
Az = N(:,4).';

%Pasamos a segundos los ms
t = t * 10^- 3;
%Hacemos que el eje empiece en 0
t = t - t(1);

%Obtenemos el periodo de muestra
Ts = t(2) - t(1);
%Obtenemos la frecuencia de muestreo
Fs = 1/Ts;
%Dibujamos cada aceleración por separado
figure(1);
subplot(3,1,1);
plot(t,Ax);
axis([t(1) t(length(t) - 1) min(Ax) max(Ax)])

subplot(3,1,2);
plot(t,Ay);
axis([t(1) t(length(t) - 1) min(Ay) max(Ay)])

subplot(3,1,3);
plot(t,Az);
axis([t(1) t(length(t) - 1) min(Az) max(Az)])


