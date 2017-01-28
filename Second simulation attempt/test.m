clc
close all
clear all


clc
close all
clear all


snr = -15:35;
M = 2;
scheme = 'BPSK';
mod = comm.PSKModulator(M,pi/M);
demod = comm.PSKDemodulator(M,pi/M);
for n = 1:length(snr)
    bits = randi([0 1],6*2^15,1);
    %channelInput = step(mod, bits);
    %bitspar = ser2par(bits,log2(M));
    channelInput = symbolGen(bits,scheme);
    
    
    hRayleighChan = comm.RayleighChannel(...
        'SampleRate',          10e3,...
        'PathDelays',          [0 1.5e-4],...
        'AveragePathGains',    [2 3],...
        'NormalizePathGains',  true,...
        'MaximumDopplerShift', 30,...
        'DopplerSpectrum',     {doppler('Gaussian',0.6), doppler('Flat')},...
        'RandomStream',        'mt19937ar with seed',...
        'Seed',                22,...
        'PathGainsOutputPort', true);
 
    % Filter the modulated data using hRayleighChan
    [chanOut1, pathGains1] = step(hRayleighChan, channelInput);
  %scatter(real(chanOut1),imag(chanOut1));
    % Now use global stream for random number generation of hRayleighChan
    release(hRayleighChan);
    hRayleighChan.RandomStream = 'Global stream';
  
    % Log the current global stream
    loggedStream = RandStream.getGlobalStream; 
  
    % Create a mt19937ar stream and designate it as current global stream
    s = RandStream('mt19937ar','Seed',22);
    RandStream.setGlobalStream(s);
  
    % Filter the modulated data using hRayleighChan for the second time
    [chanOut2, pathGains2] = step(hRayleighChan, channelInput);
    
    % Restore the logged global stream
    RandStream.setGlobalStream(loggedStream);
    y = awgn(chanOut1,snr(n));
    %channelOut = step(demod,y);
    y = y./fft(pathGains1,1,2);
    channelOut = symbolDegen(y,scheme);
    ber(n) = sum(xor(channelOut,bits))/length(bits);
end

semilogy(snr,ber,'k*')
hold on
% semilogy(snr, berawgn(snr-10*log10(1),'psk',2,'nondiff'),'LineWidth',0.5);
% semilogy(snr, berawgn(snr-10*log10(2),'psk',4,'nondiff'),'LineWidth',0.5);
% semilogy(snr, berawgn(snr-10*log10(3),'psk',8,'nondiff'),'LineWidth',0.5);
% semilogy(snr, berawgn(snr-10*log10(4),'psk',16,'nondiff'),'LineWidth',0.5);
grid;
% legend('Simulated','BPSK','QPSK','8PSK','16PSK');
return
%%
sym = [1 - 1i; -1 + 1i; -1 - 1i; 1 + 1i];
fs = 10;
n = 200;
t = 0:1/fs:n/fs;
cp = 15;

for nn = 1:length(sym)
x(nn,:) = sym(nn)*exp(1j*2*pi*fs/n*t);
plot(real(x(nn,:))+imag(x(nn,:)))
hold on

xcp(nn,:) = [x(nn,end-cp-1:end) x(nn,:)];
h = [1+1i 1+0.1i];

ycp(nn,:) = conv(xcp(nn,:),h);



hF = 2*fftshift(fft(h,1));
y(nn,:) = ycp(nn,cp:n+cp);
plot((real(y(nn,:)*exp(j*pi*0.)/hF)+imag(y(nn,:)*exp(j*pi*0.)/hF)))
yHat(nn) = mean(y(nn,:).*exp(-j*2*pi*fs/n*t));
yHat(nn) = yHat(nn)*exp(j*pi*0.4)/hF;

end
