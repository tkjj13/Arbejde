function [RecSym] = Simulation_v2_0(system,channel,symbols)


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


clc;
close all;
clear all;

clc;
close all;
clear all;

M = 4;                 % Modulation alphabet
k = log2(M);           % Bits/symbol
numSC = 128;           % Number of OFDM subcarriers
cpLen = 32;            % OFDM cyclic prefix length
maxBitErrors = 100;    % Maximum number of bit errors
maxNumBits = 1e7;      % Maximum number of bits transmitted


qpskMod = comm.QPSKModulator('BitInput',true);
qpskDemod = comm.QPSKDemodulator('BitOutput',true);

ofdmMod = comm.OFDMModulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);
ofdmDemod = comm.OFDMDemodulator('FFTLength',numSC,'CyclicPrefixLength',cpLen);

channel = comm.AWGNChannel('NoiseMethod','Variance', ...
    'VarianceSource','Input port');

errorRate = comm.ErrorRate('ResetInputPort',true);

ofdmDims = info(ofdmMod)

numDC = ofdmDims.DataInputSize(1)

frameSize = [k*numDC 1];

EbNoVec = (0:10)';
snrVec = EbNoVec + 10*log10(k) + 10*log10(numDC/numSC);

berVec = zeros(length(EbNoVec),3);
errorStats = zeros(1,3);

for m = 1:length(EbNoVec)
    snr = snrVec(m);

    while errorStats(2) <= maxBitErrors && errorStats(3) <= maxNumBits
        dataIn = randi([0,1],frameSize);                   % Generate binary data
        qpskTx = step(qpskMod,dataIn);                     % Apply QPSK modulation
        txSig = step(ofdmMod,qpskTx);                      % Apply OFDM modulation
        powerDB = 10*log10(var(txSig));                    % Calculate Tx signal power
        noiseVar = 10.^(0.1*(powerDB-snr));                % Calculate the noise variance
        rxSig = step(channel,txSig,noiseVar);              % Pass the signal through a noisy channel
        qpskRx = step(ofdmDemod,rxSig);                    % Apply OFDM demodulation
        dataOut = step(qpskDemod,qpskRx);                  % Apply QPSK demodulation
        errorStats = step(errorRate,dataIn,dataOut,0);     % Collect error statistics
    end

    berVec(m,:) = errorStats;                         % Save BER data
    errorStats = step(errorRate,dataIn,dataOut,1);         % Reset the error rate calculator
end

berTheory = berawgn(EbNoVec,'psk',M,'nondiff');

figure
semilogy(EbNoVec,berVec(:,1),'*')
hold on
semilogy(EbNoVec,berTheory)
legend('Simulation','Theory','Location','Best')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
grid on
hold off






