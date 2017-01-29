function [RecSym] = Simulation_v1_4(system,channel,symbols)

clc;
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


%%%%%%%%%%%%%%%%%
%%% VARIABLES %%%
%%%%%%%%%%%%%%%%%
nSym        = ceil(length(symbols)/system.nDSC);   % number of OFDM symbols
T           = system.fDev(system.DSCindex)^(-1);
OSF         = 5; % over samplings factor
fs          = sum(system.fDev)*OSF;
nCP         = ceil(system.CP_dur*fs*10^(-6)); % number of samples corrosponding to duration
EbN0        = channel.SNR - 10*log10(system.BPS);


%%%%%%%%%%%%%%%%%
%%% SIMULATOR %%%
%%%%%%%%%%%%%%%%%

% Generate interference symbols
missingSym = nSym*(sum(system.fDev)*T-system.nDSC)
intbits = randi([0 1],2*missingSym,1);
intSym = symbolGen(intbits,'QPSK');

% S/P
%symbols = ones(1,system.nDSC);
symbolsPar = ser2par(symbols,system.nDSC);
if (system.nFFT-system.nDSC == 0)
    intSymPar = intSym;
else
    intSymPar = ser2par(intSym,system.nFFT-system.nDSC);
end
% assign data subcarriers to OFDM subcarriers

xF = [intSymPar(:,1:system.DSCindex-1) symbolsPar intSymPar(:,system.DSCindex:end)] ;
%xF = [zeros(nSym,6) ipMod(:,[1:system.nDSC/2]) zeros(nSym,1) ipMod(:,[system.nDSC+1:system.nDSC]) zeros(nSym,5)] ;

% IFFT
sz_xF = size(xF);
xt = zeros(sz_xF(1),T*fs);%system.nFFT*OSF);
for OFDM_symbol = 1:sz_xF(1)
    for bin = 1:system.nFFT
        t = 0:1/fs:(system.fDev(system.DSCindex)^(-1))-1/fs;
        Temp = xF(OFDM_symbol,bin)*exp(1i*2*pi*(sum(system.fDev(1:bin))-system.fDev(1))*t);
        %Temp = symbolsPar(sample,bin)*exp(1i*2*pi*(system.fc(bin))*t);
        xt(OFDM_symbol,1:length(Temp)) = xt(OFDM_symbol,1:length(Temp)) + Temp;
    end
end

% Normalize ?
xt = (system.nFFT/sqrt(system.nDSC))*xt;
%xt = (system.nFFT/sqrt(system.nDSC))*ifft(fftshift(xF.')).';
% insert cyclic prefix
xt = [xt(:,end-nCP+1:end) xt];

% concatenate symbols to form long vector
xvec = reshape(xt.',1,nSym*(T*fs+nCP));

% Channel
switch (channel.type)
    case 0 % No channel
        nTap = 1;
        xht = xvec;
        hF = ones(nSym,system.nFFT);
    case 1 % multipath rayleigh fading
        nTap = ceil(channel.dt*fs);
        ht = 1/sqrt(2)*1/sqrt(nTap)*(randn(nSym,nTap) + j*randn(nSym,nTap));
        
        % computing and storing the frequency response of the channel, for use at recevier
        hF = fft(ht,system.nFFT,2);
        %hF = fftshift(fft(ht,system.nFFT,2));
        
        % convolution of each symbol with the random channel
        for jj = 1:nSym
            xht(jj,:) = conv(ht(jj,:),xvec(jj,:));
        end
    otherwise
        nTap = 1;
        xht = xvec;
        hF = ones(nSym,system.nFFT);
        %disp('No channel used');
end



yvec = awgn(xht,channel.SNR - 10*log10((T+system.CP_dur*10^(-6))/T)-10*log10(sum(system.fDev)/(system.nFFT*system.fDev(system.DSCindex))*OSF/system.BPS),'measured');

% Receiver
yht = reshape(yvec.',T*fs+nCP+nTap-1,nSym).'; % formatting the received vector into symbols
yt = yht(:,[nCP+1:nCP+T*fs]); % removing cyclic prefix



sz_yt = size(yt);
recSymPar = zeros(sz_yt(1),system.nFFT);
for sample = 1:sz_yt(1)
    for bin = 1:system.nFFT
        t = 0:1/fs:(system.fDev(bin)^(-1))-1/fs;
        Temp = exp(-1i*2*pi*(sum(system.fDev(1:bin))-system.fDev(1))*t);
        recSymPar(sample,bin) = mean(yt(sample,1:length(Temp)).*Temp);
    end
end
%recSymPar = (sqrt(system.nDSC)/system.nFFT)*fftshift(fft(yt.')).';
% equalization by the known channel frequency response
recSymPar = recSymPar./hF;

RecSym = reshape(transpose(recSymPar(:,system.DSCindex:system.DSCindex+system.nDSC-1)),length(symbols),1);
%RecSym = zeros(1,length(symbols));




