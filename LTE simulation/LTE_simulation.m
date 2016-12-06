clc
close all
clear all

rmc = struct;                      % RMC config structure
rmc.RC ='R.5';                     % Base configuration on RMC R.5
rmc.DuplexMode = 'TDD';            % User Time Division Duplex (TDD)
rmc.TotSubframes = 1;              % Configure a single subframe

% Generate the base configuration from RMC R.5 and amend to set the
% parameters required for the Port7-14 transmission scheme. Note that if
% the standard defined RMCs using Port7-14 transmission scheme supported by
% lteRMCDL is used, these parameters will be pre-configured.
rmc = lteRMCDL(rmc);
rmc.NDLRB = 25;                    % 25 Resource Blocks
rmc.NCellID = 10;                  % Cell identity 10
rmc.PDSCH.TxScheme = 'Port7-14';   % Up to 8 layer transmission, ports 7-14
rmc.PDSCH.NLayers = 2;             % 2 transmission layers for beamforming
rmc.PDSCH.NSCID = 0;               % Scrambling identity 0
rmc.CSIRefP = 8;                   % 8 CSI-RS ports
rmc.CSIRSConfig = 0;               % CSI-RS configuration 0
rmc.CSIRSPeriod = 'On';            % Configure CSI-RS always 'on'
rmc.ZeroPowerCSIRSPeriod = 'Off';  % Configure Zero Power CSI-RS 'off'
rmc.PDSCH.PRBSet = (4:8).';        % 5 allocated RBs
rmc.PDSCH.PMIMode = 'Wideband';    % Wideband PMI mode
rmc.PDSCH.CSI = 'On';              % CSI scaling of soft bits
% Codebook subset definition allowing all codebook entries
rmc.PDSCH.CodebookSubset = '0x1FFFFFFFFFFFFFFFFFFFFFFFFFFF';




channel = struct;                   % Channel config structure
channel.Seed = 8;                   % Random channel seed
channel.NRxAnts = 3;                % 3 receive antennas
channel.DelayProfile = 'EVA';       % Delay profile
channel.DopplerFreq = 5.0;          % Doppler frequency in Hz
channel.MIMOCorrelation = 'Medium'; % Multi-antenna correlation
channel.NTerms = 16;                % Oscillators used in fading model
channel.ModelType = 'GMEDS';        % Rayleigh fading model type
channel.InitTime = 0.0;             % Initial time
channel.InitPhase = 'Random';       % Random initial phases
channel.NormalizePathGains = 'On';  % Normalize delay profile power
channel.NormalizeTxAnts = 'On';     % Normalize for transmit antennas



cec = struct;                       % Channel estimation config structure
cec.PilotAverage = 'UserDefined';   % Type of pilot symbol averaging
cec.FreqWindow = 1;                 % Frequency window size (special mode)
cec.TimeWindow = 2;                 % Time window size (special mode)
cec.InterpType = 'Cubic';           % 2D interpolation type
cec.InterpWindow ='Centered';       % Interpolation window type
cec.InterpWinSize = 1;              % Interpolation window size






% Transmit without then with CSI-RS-based PMI feedback
for csirsFeedback = 0:1

    % Configure random number generators
    rng('default');

    % Configure PDSCH substructure with transmission beamforming matrix W.
    % In the first iteration of the loop transmit each layer on one of the
    % 8 antennas. In the second iteration, transmit the layers on 2 beams
    % matched to the channel response using CSI-RS-based PMI feedback. The
    % PMI fed back to the second iteration is calculated at the end of the
    % first
    if ~csirsFeedback
        rmc.PDSCH.W = [1 0 0 0 0 0 0 0; ...
                       0 0 0 0 1 0 0 0]/sqrt(2);
    else
        rmc.PDSCH.W = lteCSICodebook(rmc.PDSCH.NLayers, ...
                                     rmc.CSIRefP, [PMI(1) PMI(2)]).';
    end

    % Generate transmission with PDSCH with the beamforming matrix W, onto
    % 1st of 8 antenna planes (note that CellRefP = 1 for this RMC). The
    % transmitted grid contains UE-specific reference signal (UE-RS / DMRS)
    % for channel estimation and CSI-RS reference signal for PMI selection
    [~, txGrid, rmcinfo] = lteRMCDLTool(rmc, [1;0;0;1]);
    channel.SamplingRate = rmcinfo.SamplingRate;

    % OFDM modulation. The additional 25 samples added to the end of the
    % waveform are to cover the range of delays expected from the channel
    % modeling (a combination of implementation delay and channel delay
    % spread)
    [txWaveform, ofdmDims] = lteOFDMModulate(rmc, txGrid, 0);
    txWaveform = [txWaveform; zeros(25, size(txWaveform,2))]; %#ok

    % Fading channel
    rxWaveform = lteFadingChannel(channel, txWaveform);

    % Create and apply additive white Gaussian noise
    if ~csirsFeedback
        SNRdB = 27;
        SNR = 10^(SNRdB/20);
        N = 1/(sqrt(2.0*rmc.CSIRefP*double(ofdmDims.Nfft))*SNR);
        v = N*complex(randn(size(rxWaveform)), randn(size(rxWaveform)));
    end
    rxWaveform = rxWaveform + v;

    % Perform synchronization
    offset = lteDLFrameOffset(rmc, rxWaveform);
    rxWaveform = rxWaveform(1+offset:end, :);

    % Perform OFDM demodulation on the received data to recreate the
    % resource grid
    rxGrid = lteOFDMDemodulate(rmc, rxWaveform);

    % Channel estimation using the UE-specific DMRS for PDSCH reception
    cec.Reference = 'DMRS';
    [hest, nest] = lteDLChannelEstimate(rmc, rmc.PDSCH, cec, rxGrid);

    % Equalize (back to layers) and demodulate the PDSCH.
    % Extract REs corresponding to the 2 layers of the PDSCH from the given
    % subframe across all receive antennas and channel estimates.
    ind = ltePDSCHIndices(rmc, rmc.PDSCH, rmc.PDSCH.PRBSet);
    [pdschRx, pdschHest] = lteExtractResources(ind, rxGrid, hest);
    [rxBits, rxSymbols] = ltePDSCHDecode(rmc, rmc.PDSCH, ...
        pdschRx, pdschHest, nest);

    % Compute singular values of the channel and calculate SNR
    H = squeeze(mean(pdschHest));
    d = svd(H);

    % Print singular values and effective channel SNR
    if csirsFeedback
        label = '8 antenna transmission with CSI-RS-based PMI feedback';
    else
        label = '8 antenna transmission, 1 antenna for each layer';
    end
    fprintf('%s:\n\n', label);
    svdb = sprintf('     %0.2fdB', 20*log10(d));
    fprintf('       Channel singular values:%s\n', svdb);
    fprintf('         Effective channel SNR:     %0.2fdB\n', ...
        SNRdB+10*log10(rmc.PDSCH.NLayers)+10*log10(sum(d.^2)));

    % Regenerate PDSCH from hard bit decisions and demodulate to estimate
    % transmitted symbols
    remod = ltePDSCH(rmc, rmc.PDSCH, rxBits{1}>0);
    [rxBitsRef, rxSymbolsRef] = ltePDSCHDecode(rmc, rmc.PDSCH, remod);

    % Use EVM measurement to estimate SNR
    EVM = comm.EVM;
    evmRMS = step(EVM,rxSymbols{1}, rxSymbolsRef{1});
    SNRest = 20*log10(1/(evmRMS/100));
    fprintf('SNR estimate from receiver EVM:     %0.2fdB\n\n',SNRest);

    % Now compute PMI (via CSI-RS) for use in second iteration. Channel
    % realization remains the same
    if ~csirsFeedback
        % Channel estimation via CSI-RS for PMI selection
        cec.Reference = 'CSIRS';
        [hestPMI, nestPMI] = lteDLChannelEstimate(rmc, rmc.PDSCH, ...
            cec, rxGrid);
        % PMI selection
        PMI = ltePMISelect(rmc, rmc.PDSCH, hestPMI, nestPMI);
    end

    % Plot received constellation
    figure(csirsFeedback+1);
    plot(rxSymbols{1}, 'o', 'MarkerEdgeColor', [0.75 0 0], ...
        'MarkerFaceColor', [1 0.25 0.25], 'MarkerSize',3);
    axis([-1.25 1.25 -1.25 1.25]);
    title(label);

end