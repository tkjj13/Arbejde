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
number_of_symbols = 200;
T = 1E-3;
N = 1024;
CP_rate = 0;
EbN0 = 15;
omega = 1;
sheme = 'QPSK';
BPS = 2;    % bits per symbol

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
Bits = randi([0 1],BPS*number_of_symbols*number_of_cariers,1);
%Bits = repmat([1 1 0 1 1 0 0 0]',10,1);
[symbols] = symbolGen(Bits, sheme);

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
h = sqrt(-log(rand(1,length(xSer)))*omega);
y1 = h.*xSer;     % rayleigh fading



% noise

%/***	snr=Eb_No-10log(BW/rb)	***/	
snr_mark = EbN0 - 10*log10(N);	
%/***	=eff./10^(snr_mark/10)		***/	
sigma_i_mark = 1/10^(snr_mark/10);	
%/***	Add noise	***/	
x_xi = rand(1,length(y1));			% /***	random	number	between	0	and	1	***/	
x_psi = rand(1,length(y1));	
% if(x_xi>=1.0)
%     x_xi = 0.99999
%     disp('Rand	overflow!!!\n');
% end
xi = sqrt(-2*sigma_i_mark*log10(1.0-x_xi));	
y = y1+(xi.*cos(2*pi*x_psi)+1i*xi.*sin(2*pi*x_psi));	%/*ri(t)=si(t)+ni(t)*/	 %/*rq(t)=sq(t)+nq(t)*/	
	



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
recSym = reshape(transpose(recSymPar),sample*bin,1);
scatter(real(recSym),imag(recSym));

recBits = symbolDegen(recSym,sheme);

Err =  sum(xor(recBits,Bits));
BER = (Err/length(Bits))

