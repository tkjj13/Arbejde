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
OSF         = 5; % over samplings factor
for n = 1:length(system.fDev)
    fc(n)   = sum(system.fDev(1:n));
end
fs          = sum(system.fDev)*OSF;
nCP         = ceil(system.CP_dur*fs*10^(-6)); % number of samples corrosponding to duration
EbN0        = channel.SNR - 10*log10(system.BPS);
fDevMax     = max(system.fDev);
TMin        = fDevMax^(-1);
fDevDSC     = system.fDev(system.DSCindex);
TDSC        = fDevDSC^(-1);
%%%%%%%%%%%%%%%%%
%%% SIMULATOR %%%
%%%%%%%%%%%%%%%%%

% Generate interference symbols
missingSym = nSym*(sum(system.fDev)*TDSC-system.nDSC);
intbits = randi([0 1],2*missingSym,1);
intSym = symbolGen(intbits,'QPSK');

% S/P
%symbols = ones(1,system.nDSC);
symbolsPar = ser2par(symbols,system.nDSC);

% assign data subcarriers to OFDM subcarriers

xF = zeros(fDevMax/fDevDSC*nSym,system.nFFT);
index = 0;
for n = 1:system.nFFT
    if ((n < system.DSCindex) || (n > system.DSCindex+system.nDSC-1))
        xF(:,n) = reshape(transpose(repmat(intSym(index+1:index+nSym*system.fDev(n)/fDevDSC),1,fDevMax/system.fDev(n))),fDevMax/fDevDSC*nSym,1);
        index = index+nSym*system.fDev(n)/fDevDSC;
    else
        xF(:,n) = reshape(transpose(repmat(symbolsPar(:,n-system.DSCindex+1),1,fDevMax/system.fDev(n))),fDevMax/fDevDSC*nSym,1);
    end
end
    %xF = [reshape(intSymPar(:,1:system.DSCindex-1) symbolsPar intSymPar(:,system.DSCindex:end)] ;
    %xF = [zeros(nSym,6) ipMod(:,[1:system.nDSC/2]) zeros(nSym,1) ipMod(:,[system.nDSC+1:system.nDSC]) zeros(nSym,5)] ;

% IFFT
sz_xF = size(xF);
xt = zeros(sz_xF(1), sz_xF(2), TMin*fs);

for FDM_symbol = 1:sz_xF(1)
    for bin = 1:system.nFFT
       t = (FDM_symbol-1)*TMin:1/fs:FDM_symbol*TMin-1/fs;
%         t = 0:1/fs:system.fDev(1)^(-1)-1/fs;
%         xt(FDM_symbol,:,:) = diag(xF(FDM_symbol,:))*exp(1i*2*pi*fc'*t);
        xt(FDM_symbol,bin,:) = xF(FDM_symbol,bin)*exp(1i*2*pi*fc(bin)*t);      
%         Temp = symbolsPar(sample,bin)*exp(1i*2*pi*(system.fc(bin))*t);
%         Temp2 = xF(FDM_symbol,bin)*exp(1i*2*pi*(sum(system.fDev(1:bin))-system.fDev(1))*t);
%         xt2(FDM_symbol,1:length(Temp2)) = xt2(FDM_symbol,1:length(Temp2)) + Temp2;
    end
end
xt2(:,:) = sum(xt,2); 
% Normalize ?
xt = (system.nFFT/sqrt(system.nDSC))*xt;
%xt = (system.nFFT/sqrt(system.nDSC))*ifft(fftshift(xF.')).';

% insert cyclic prefix and concatenate
for bin = 1:system.nFFT
    Temp = reshape(xt(:,bin,:),sz_xF(1)*system.fDev(bin)/fDevMax,fDevMax/system.fDev(bin)*TMin*fs);
    Temp = [Temp(:,end-ceil(nCP*fDevDSC/system.fDev(bin))+1:end) Temp];
    sz_temp = size(Temp);
    xCP(bin,:) = reshape(transpose(Temp),1,sz_temp(1)*sz_temp(2));
end
xvec = sum(xCP,1);
%xt2 = [xt2(:,end-nCP+1:end) xt2];

% concatenate symbols to form long vector
%xvec = reshape(xt2.',1,nSym*(TDSC*fs+nCP));

% Channel
switch (channel.type)
    case 0 % No channel
        nTap = 1;
        xht = xvec;
        hF = ones(1,system.nFFT);
    case 1 % multipath rayleigh fading
        chan = rayleighchan(1/fs,100);
        xht = filter(chan,xvec);
        
    otherwise
        nTap = 1;
        xht = xvec;
        hF = ones(1,system.nFFT);
        %disp('No channel used');
end



yvec = awgn(xht,channel.SNR,'measured');
% yvec = xht;
% image rejection
% Wn = sum(system.fDev(1:system.nDSC-5))/fs;
% [b, a] = butter(30, Wn);
% yvec = filter(b,a, yvec);

% Receiver
yht = reshape(transpose(yvec(nTap:end)),TDSC*fs+nCP,nSym).'; % formatting the received vector into symbols
yt = yht(:,[nCP+1:nCP+TDSC*fs]); % removing cyclic prefix



sz_yt = size(yt);
recSymPar = zeros(sz_yt(1),system.nFFT);
for sample = 1:sz_yt(1)
    for bin = 1:system.nFFT
        t = 0:1/fs:system.fDev(bin)^(-1)-1/fs;
        Temp = exp(-1i*2*pi*(sum(system.fDev(1:bin)))*t);
        recSymPar(sample,bin) = mean(yt(sample,1:length(Temp)).*Temp);
    end
end
%recSymPar = (sqrt(system.nDSC)/system.nFFT)*fftshift(fft(yt.')).';
% equalization by the known channel frequency response
%recSymPar = recSymPar./repmat(hF,100,1);

RecSym = reshape(transpose(recSymPar(:,system.DSCindex:system.DSCindex+system.nDSC-1)),length(symbols),1);
%RecSym = zeros(1,length(symbols));
