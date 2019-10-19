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

%Obtenemos la amgnitud del vector 3D
tam_array = length(t); 
A = zeros(1,tam_array);

for i = 1:tam_array -1
    
    %x = Ax(i) - Ax(i + 1);
    %y = Ay(i) - Ay(i + 1);
    %z = Az(i) - Az(i + 1);
    %A(1,i) = sqrt(x^2 + y^2 + z^2);
    
    A(1,i) = sqrt(Ax(1,i)^2 + Ay(1,i)^2 + Az(1,i)^2);
end


%%Filtramos
%%Definimos filtro IIR PASO-ALTO
f0 = 1;
[b a] = butter(7, 2*f0/Fs, 'high');

%%Definimos filtro PASO-BAJO
f1 = 20;
[b1 a1] = butter(4, 2*f1/Fs, 'low');

%%Filtramos

A_filtrada = filtfilt(b,a,A);
A_filtrada = filtfilt(b1, a1, A_filtrada);

%Integral TRAPEZOIDAL
V = zeros(1,tam_array);

for i = 1:tam_array - 1
    V(1,i) = trapz(A_filtrada(1,i:i+1));
end

P = zeros(1,tam_array);

for i = 1:tam_array - 1
    P(1,i) = trapz(V(1,i:i+1));
end

%Respuesta en frecuencia
A1 = A.*hann(length(A))';
L = nextpow2(length(A1));
[X1, f] = freqz(A1, 1, 2^L, Fs);
P1 = 1/2^L*abs(X1).*2;
 
A2_filtrada = A_filtrada.*hann(length(A_filtrada))';
L = nextpow2(length(A2_filtrada));
[X2, f] = freqz(A2_filtrada, 1,2^L,Fs);
P2 = 1/2^L*abs(X2).*2;

%Dibujamos dominio de la frecuencia
figure(1);
subplot(2,1,1);
plot(f,P1);

subplot(2,1,2);
plot(f,P2);

%Dibujamos cada aceleración por separado
figure(2)
subplot(6,1,1);
plot(t,Ax);
axis([t(1) t(length(t) - 1) min(Ax) max(Ax)])

subplot(6,1,2);
plot(t,Ay);
axis([t(1) t(length(t) - 1) min(Ay) max(Ay)])

subplot(6,1,3);
plot(t,Az);
axis([t(1) t(length(t) - 1) min(Az) max(Az)])

subplot(6,1,4);
plot(t,A_filtrada);
hold on;
plot(t,A);

subplot(6,1,5)
plot(t,V);

subplot(6,1,6)
plot(t,P);
%axis([t(1) t(length(t) - 1) min(A_filtrada) max(A_filtrada)])


