% clear
% clc
fprintf('Start of script %s\n', datestr(now,'HH:MM:SS.FFF'))
params
%% Create Transmitter and Receiver objects
transmitter = comm.SDRTxPluto('RadioID', 'ip:192.168.3.1');

%% Sampling
transmitter.CenterFrequency = Fc;
transmitter.BasebandSampleRate = Fs;

%% Transmit Signal
signalIQ_pre = LoRa_Tx(message,BW,SF,Power,Fs,Fc - fc) ;
%No_symb = length(signalIQ_pre);
%n = 1/sqrt(2)*[randn(1,No_symb) + 1i*randn(1,No_symb)];
%h = generate_patzold(kappa_0, sigma_0, s_param, m, kappa_c, N1, N2, f_max, No_symb);
%signalIQ = h'.*signalIQ_pre + 10^(-(SNR)/20)*n';
signalIQ = signalIQ_pre;
fprintf('Right before transmitting %s\n', datestr(now,'HH:MM:SS.FFF'))
pause(0.1)
transmitter(signalIQ);
fprintf('Right after transmitting %s\n', datestr(now,'HH:MM:SS.FFF'))

Sxx = 10*log10(rms(signalIQ).^2) ;
disp(['Transmit Power   = ' num2str(Sxx) ' dBm'])
%% Plots
figure()
spectrogram(signalIQ_pre,500,0,500,Fs,'yaxis','centered')
figure()
obw(signalIQ,Fs) ;

