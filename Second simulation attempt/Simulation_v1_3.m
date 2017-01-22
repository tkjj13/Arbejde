function [RecSym] = Simulation_v1_3(system,channel,symbols)


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



% system.fDev = [repmat(1000,1,5) repmat(2000,1,5)];
system.fc = [0 10000 20000 30000 40000 60000 80000 100000 120000 140000];
% system.N = 1024;
 system.CP_rate = 0.75;


%%%%%% Calculating secondary parameters

T = system.fDev.^(-1);
Tmax = max(T);
number_of_cariers = length(system.fDev);
fs = system.N/Tmax;

number_of_symbols = length(symbols)/length(system.fDev);
%%% USED ONLY FOR GENERATING SYMBOLS %%%


% S/P
symbolsPar = ser2par(symbols,number_of_cariers);

% IFFT
%t = linspace(0,Tmax,N);             % consider if it should be t = linspace(0,Tmax-Tmax/N,N);
sz_par = size(symbolsPar);
x = zeros(sz_par(1),system.N);
for sample = 1:sz_par(1)
    for bin = 1:number_of_cariers
        t = 0:1/fs:T(bin)-1/fs;
        %Temp = symbolsPar(sample,bin)*exp(1i*2*pi*(sum(fDev(1:bin))-fDev(1))*t);
        Temp = symbolsPar(sample,bin)*exp(1i*2*pi*(system.fc(bin))*t);
        x(sample,1:length(Temp)) = x(sample,1:length(Temp)) + Temp;
    end
end


% x = N*ifft(symbolsPar,N,2);


% Insert CP
sz_x = size(x);

samples = system.CP_dur*0.001*fs;
reps = floor(samples/system.N);
rest = floor(mod(samples,system.N));

xCp = repmat(x,1,reps+1);
xCp = [xCp(:,end-rest:end) xCp(:,2:end)];
%xCp = [x(:,round(sz_x(2)*(1-system.CP_rate)):end) x(:,2:end)];

% P/S
sz_xCp = size(xCp);
xSer = reshape(transpose(xCp),1,sz_xCp(1)*sz_xCp(2));


% Channel
switch (channel.type)
    case 0 % singlepath rayleigh fading
        Omega = 1;
        h = sqrt(-log(rand(1,length(xSer)))*1);
        y1 = h.*xSer;     % rayleigh fading
    
    case 1 % multipath rayleigh fading
        offset = ceil(channel.dt*fs);
        y1 = zeros(1,length(xSer)+(length(channel.h)-1)*offset);
    
        for k = 0:length(channel.h)-1
            s = zeros(size(y1));
            s(k*offset+1:end-(length(channel.h)-1-k)*offset) = channel.h(k+1)*xSer;
            y1 = y1+s;
        end
    otherwise 
        y1 = xSer;
        disp('No channel used');
end


% noise

%	snr=Eb_No-10log(BW/rb)		
snr_mark = channel.EbN0 - 10*log10(system.N);	
%	=eff./10^(snr_mark/10)			
sigma_i_mark = 1/10^(snr_mark/10);	
%	Add noise		
x_xi = rand(1,length(y1));			%	random	number	between	0	and	1
x_psi = rand(1,length(y1));	

xi = sqrt(-2*sigma_i_mark*(x_xi));	
y = y1+(xi.*cos(2*pi*x_psi)+1i*xi.*sin(2*pi*x_psi));	%ri(t)=si(t)+ni(t)	 rq(t)=sq(t)+nq(t)	
	



% % image rejection
% Wn = 0.6;
% [b, a] = butter(10, Wn);
% y = filter(b,a, ySer);



% Correct CFO
% S/P
yParCP = transpose(reshape(y(1:length(xSer)),system.N+samples,number_of_symbols));

% Remove CP
yPar = yParCP(:,end-system.N:end);



% FFT

sz_yPar = size(yPar);
recSymPar = zeros(sz_yPar(1),number_of_cariers);
for sample = 1:sz_yPar(1)
    for bin = 1:number_of_cariers
        t = 0:1/fs:T(bin)-1/fs;
        %Temp = exp(-1i*2*pi*(sum(fDev(1:bin))-fDev(1))*t);
        Temp = exp(-1i*2*pi*(system.fc(bin))*t);
%         if bin ~= 1
%             Temp2 = yPar(sample,1:length(Temp));
%             for k = 1:bin-1
%                 t = 0:1/fs:T(bin)-1/fs;
%                 %Temp3 = recSymPar(sample,bin)*exp(1i*2*pi*(sum(fDev(1:bin))-fDev(1))*t);
%                 Temp3 = recSymPar(sample,bin)*exp(1i*2*pi*system.fc(k)*t);
%                 Temp2 = Temp2 - Temp3;
%             end
%             recSymPar(sample,bin) = mean(Temp2.*Temp);
%         else
%             recSymPar(sample,bin) = mean(yPar(sample,1:length(Temp)).*Temp);
%         end
       recSymPar(sample,bin) = mean(yPar(sample,1:length(Temp)).*Temp);
        RecSym = reshape(transpose(recSymPar),sz_yPar(1)*number_of_cariers,1);
%plot(t,real(yPar(sample,1:length(Temp)).*Temp),t,imag(yPar(sample,1:length(Temp)).*Temp));
    end
end

 


% Phase track
% P/S
RecSym = reshape(transpose(recSymPar),sample*bin,1);
%scatter(real(RecSym),imag(RecSym));



