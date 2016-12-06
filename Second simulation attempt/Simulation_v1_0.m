clc;
close all;
clear all;



% input
% S/P
% IFFT
% Insert CP
% P/S
% D/A
% Channel
% noise
% A/D
% Correct CFO
% S/P
% Remove CP
% FFT
% Phase track
% P/S



number_of_cariers = 10;
T = 1E-3;
N = 1000;
CP_rate = 0;
fc = 1E6;
H = [1];
H_dt = T/2;


% input
Bits = randi([0 1],10*number_of_cariers,1);
[symbols] = symbolGen(Bits, 'QPSK');

% S/P
symbolsPar = ser2par(symbols,number_of_cariers);

% IFFT
t = T/N:T/N:T;
Length_t = length(t);
sz_par = size(symbolsPar);
x = zeros(sz_par(1),Length_t);
for sample = 1:sz_par(1)
    for bin = 1:number_of_cariers
        x(sample,:) = x(sample,:) + symbolsPar(sample,bin)*exp(1i*2*pi*(bin-1)*t/T);
    end
end


% Insert CP
xCp = [x(:,round(Length_t*(1-CP_rate)):end) x(:,2:end)];

% P/S
sz_xCp = size(xCp);
xSer = reshape(transpose(xCp),1,sz_xCp(1)*sz_xCp(2));

% scale to RF  
t = T/N:T/N:T*length(xSer)/N;
RF_tx = real(xSer).*sqrt(2/T).*cos(2*pi*fc*t)+imag(xSer).*sqrt(2/T).*sin(2*pi*fc*t);

% Channel
sz_H = size(H);
RF_channel = zeros(1,length(RF_tx)+(sz_H(1)-1)*H_dt*N/T);
for n = 0:sz_H(1)-1;
    RF_channel = RF_channel + [zeros(1,n*ceil(H_dt*N/T)) RF_tx zeros(1,(sz_H(1)-1-n)*ceil(H_dt*N/T))]*H(n+1);
end

% noise



% scale to complex baseband
t = T/N:T/N:T*length(RF_channel)/N;
ySer = 2*RF_channel./sqrt(2/T).*exp(1i*2*pi*fc*t);

% image rejection
Wn = 0.6;
[b, a] = butter(10, Wn);
y = filter(b,a, ySer);


figure;
plot(real(xSer))
hold on
plot(real(y))

figure;
plot(imag(xSer))
hold on
plot(imag(y))

figure
plot(abs(fft(T*real(y))));


% Correct CFO
% S/P
%yParCP = ser2par(transpose(y),N*(1+CP_rate));
yParCP = transpose(reshape(y,N*(1+CP_rate),5));



% Remove CP
yPar = yParCP(:,round(N*CP_rate)+1:end);

% figure
% plot(real(x(1,:)))
% hold on
% plot(real(yPar(1,:)))
% return



% FFT
t = T/N:T/N:T;
Length_t = length(t);

sz_yPar = size(yPar);
recSymPar = zeros(sz_yPar(1),number_of_cariers);
for sample = 1:sz_yPar(1)
    for bin = 1:number_of_cariers
        recSymPar(sample,bin) = mean((yPar(sample,:)).*exp(-1i*2*pi*(bin-1)*t/T));
    end
end

% Phase track
% P/S
symbolsPar-recSymPar



