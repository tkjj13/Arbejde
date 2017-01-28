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
OSF         = 1; % over samplings factor
fs          = sum(system.fDev)*OSF;
nCP         = ceil(system.CP_dur*fs*10^(-6)); % number of samples corrosponding to duration
EbN0        = channel.SNR - 10*log10(system.BPS*(1-system.fDev(system.DSCindex)*system.CP_dur*10^(-6)));



%%%%%%%%%%%%%%%%%
%%% SIMULATOR %%%
%%%%%%%%%%%%%%%%%

% S/P
%symbols = ones(1,system.nDSC);
symbolsPar = ser2par(symbols,system.nDSC);

% assign data subcarriers to OFDM subcarriers
xF = [symbolsPar zeros(nSym,system.nFFT-system.nDSC)] ;
%xF = [zeros(nSym,6) ipMod(:,[1:system.nDSC/2]) zeros(nSym,1) ipMod(:,[system.nDSC+1:system.nDSC]) zeros(nSym,5)] ;

% IFFT
sz_xF = size(xF);
xt = zeros(sz_xF(1),system.nFFT*OSF);
for OFDM_symbol = 1:sz_xF(1)
    for bin = 1:system.nFFT
        t = 0:1/fs:(system.fDev(bin)^(-1))-1/fs;
        Temp = xF(OFDM_symbol,bin)*exp(1i*2*pi*(sum(system.fDev(1:bin))-system.fDev(1))*t);
        %Temp = symbolsPar(sample,bin)*exp(1i*2*pi*(system.fc(bin))*t);
        xt(OFDM_symbol,1:length(Temp)) = xt(OFDM_symbol,1:length(Temp)) + Temp;
    end
end

% Normalize ?
%xt = (system.nFFT/sqrt(system.nDSC))*xt;
xt = (system.nFFT/sqrt(system.nDSC))*ifft(fftshift(xF.')).';
% insert cyclic prefix
xt = [xt(:,end-nCP+1:end) xt];

% concatenate symbols to form long vector

% Channel
switch (channel.type)
    case 0 % No channel
        nTap = 1;
        xht = xt;
        hF = ones(nSym,system.nFFT);
    case 1 % multipath rayleigh fading
        nTap = ceil(channel.dt*fs);
        ht = 1/sqrt(2)*1/sqrt(nTap)*(randn(nSym,nTap) + j*randn(nSym,nTap));

        % computing and storing the frequency response of the channel, for use at recevier
        hF = fft(ht,system.nFFT,2);
        hF = fftshift(fft(ht,system.nFFT,2));

        % convolution of each symbol with the random channel
        for jj = 1:nSym
           xht(jj,:) = conv(ht(jj,:),xt(jj,:));
        end
    otherwise 
        nTap = 1;
        xht = xt;
        hF = ones(nSym,system.nFFT);
        %disp('No channel used');
end

% concatenate symbols to form long vector
xvec = reshape(xht.',1,nSym*(OSF*system.nFFT+nCP+nTap-1));

%    % Gaussian noise of unit variance, 0 mean
%    nt = 1/sqrt(2)*[randn(1,length(xvec)) + j*randn(1,length(xvec))];
% 
%    % Adding noise, the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
%    yvec = sqrt((system.nFFT+nCP)/system.nFFT)*xvec + 10^(EbN0/20)*nt;

snr_mark = EbN0 - 10*log10(OSF*system.nFFT);	
%	=eff./10^(snr_mark/10)			
sigma_i_mark = 1/10^(snr_mark/10);	
%	Add noise		
x_xi = rand(1,length(xvec));			%	random	number	between	0	and	1
x_psi = rand(1,length(xvec));	

xi = sqrt(-2*sigma_i_mark*(x_xi));	
yvec = xvec+(xi.*cos(2*pi*x_psi)+1i*xi.*sin(2*pi*x_psi));	%ri(t)=si(t)+ni(t)	 rq(t)=sq(t)+nq(t)
   
   
% Receiver
   yht = reshape(yvec.',OSF*system.nFFT+nCP+nTap-1,nSym).'; % formatting the received vector into symbols
   yt = yht(:,[nCP+1:nCP+system.nFFT*OSF]); % removing cyclic prefix

   
   
sz_yt = size(yt);
recSymPar = zeros(sz_yt(1),system.nFFT);
for sample = 1:sz_yt(1)
    for bin = 1:system.nFFT
        t = 0:1/fs:(system.fDev(bin)^(-1))-1/fs;
        Temp = exp(-1i*2*pi*(sum(system.fDev(1:bin))-system.fDev(1))*t);
        recSymPar(sample,bin) = mean(yt(sample,1:length(Temp)).*Temp);
    end
end   
recSymPar = (sqrt(system.nDSC)/system.nFFT)*fftshift(fft(yt.')).';
% equalization by the known channel frequency response
   recSymPar = recSymPar./hF;

RecSym = reshape(transpose(recSymPar(:,system.DSCindex:system.DSCindex+system.nDSC-1)),length(symbols),1);
%RecSym = zeros(1,length(symbols));




