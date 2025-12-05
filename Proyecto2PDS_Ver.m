% Filter Design 2026
clf;
clear all;
format long g

M = 2000;
N = 250;
fo = 15;          % frecuencia central
D = 5;             % ancho de banda del stopband
n = 2;             % orden del filtro

fo1 = fo + fo;
fo2 = fo + 2*fo;
fo3 = fo + 3*fo;
fo4 = fo + 4*fo;
fo5 = fo + 5*fo;

A = 1;

% Señal con 6 tonos
nvec = 1:M;
x = A*cos(2*pi*fo*nvec/M) + ...
    A*cos(2*pi*fo1*nvec/M)  + ...
    A*cos(2*pi*fo2*nvec/M)  + ...
    A*cos(2*pi*fo3*nvec/M)  + ...
    A*cos(2*pi*fo4*nvec/M)  + ...
    A*cos(2*pi*fo5*nvec/M);

Xe = abs(fft(x));

% Frecuencias normalizadas (Butterworth) para fo1
wL = (fo1 - D) * 2 / M;   % límite inferior normalizado
wH = (fo1 + D) * 2 / M;   % límite superior normalizado
w  = [wL wH];
% Para f03
wL3 = (fo3 - D) * 2 / M;   % límite inferior normalizado
wH3 = (fo3 + D) * 2 / M;   % límite superior normalizado
w3  = [wL3 wH3];

% Diseño del filtro stopband
[b,a] = butter(n, w, 'stop');
disp('Coeficientes del filtro 1:')
disp(b)
disp(a)
% Diseño del filtro stopband para fo3
[b3, a3] = butter(n, w3, 'stop');
disp('Coeficientes del filtro 3:')
disp(b3)
disp(a3)
% Filtrado con el primer filtro
y = filter(b, a, x);
% Filtrado con el segundo filtro
y = filter(b3, a3, y);
Ye = abs(fft(y));

 % Respuesta en frecuencia
[H, wf] = freqz(b, a, N/2);

%% Gráficas
figure(1)
subplot(4,1,1), plot(x); title('Señal original');
subplot(4,1,2), plot(Xe, 'r'); title('FFT de la señal original');
subplot(4,1,3), plot(y); title('Señal filtrada');
subplot(4,1,4), plot(Ye); title('FFT de la señal filtrada');

grid on;
