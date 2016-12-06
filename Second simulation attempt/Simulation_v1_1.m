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
number_of_symbols = 5;
T = 1E-3;
N = 1024;
CP_rate = 0;
EbN0 = 3;

%%%%%%%%%%% NOT USED %%%%%%%%%%%%%
k = 1;
H = ones(1,k*number_of_cariers);
H_df = 1/(k*T);

if length(H)*H_df ~= number_of_cariers/T
    disp('Channel characterization is not covering the whole channel');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fs = N/T;

% input
Bits = randi([0 1],2*number_of_symbols*number_of_cariers,1);
[symbols] = symbolGen(Bits, 'QPSK');

% S/P
symbolsPar = ser2par(symbols,number_of_cariers);

% IFFT
% t = 1/fs:1/fs:T;
% Length_t = length(t);
% sz_par = size(symbolsPar);
% x = zeros(sz_par(1),Length_t);
% for sample = 1:sz_par(1)
%     for bin = 1:number_of_cariers
%         x(sample,:) = x(sample,:) + symbolsPar(sample,bin)*exp(1i*2*pi*(bin-1)*t/T);
%     end
% end

x = N*ifft(symbolsPar,N,2);


% Insert CP
sz_x = size(x);
xCp = [x(:,round(sz_x(2)*(1-CP_rate)):end) x(:,2:end)];

% P/S
sz_xCp = size(xCp);
xSer = reshape(transpose(xCp),1,sz_xCp(1)*sz_xCp(2));


% Channel
h = sqrt(-log(rand(1,length(xSer))));
y = h.*xSer;     % rayleigh fading



% noise




% % image rejection
% Wn = 0.6;
% [b, a] = butter(10, Wn);
% y = filter(b,a, ySer);



% Correct CFO
% S/P
yParCP = transpose(reshape(y,N*(1+CP_rate),number_of_symbols));

% Remove CP
yPar = yParCP(:,round(N*CP_rate)+1:end);




% FFT
t = 0:1/fs:T-1/fs;
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
recSym = reshape(recSymPar,sample*bin,1);
scatter(real(recSym),imag(recSym));



