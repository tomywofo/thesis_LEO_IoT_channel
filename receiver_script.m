clear
clc

fprintf('Start of script %s\n', datestr(now,'HH:MM:SS.FFF'))
params
%% Create Transmitter and Receiver objects
receiver = comm.SDRRxPluto('RadioID', 'ip:192.168.2.1');

%% Sampling

receiver.CenterFrequency = Fc;
receiver.BasebandSampleRate = Fs;
%% Received Signal
fprintf('Right before receiving %s\n', datestr(now,'HH:MM:SS.FFF'))
rxLogNoOverflow = dsp.SignalSink;
rxLogDataValid = dsp.SignalSink;
rxLogAllData = dsp.SignalSink;
pause(0.2)
for approx_len = 1:100
    [data,datavalid,overflow] = receiver();
% Check for overflow of received samples.
    if (overflow) % dropped samples
        disp('samples dropped');
    else
        rxLogNoOverflow(data);
        disp('invalid data logged')
    end
% Check for overflow and validity of received data.
    if ~(overflow) % no dropped samples
        if ~(datavalid) % received desired data
             rxLogDataValid(data);
             disp('valid data logged')
        end
    else
        disp('no valid data received');
    end
    rxLogAllData(data);
end
fprintf('Right after receiving %s\n', datestr(now,'HH:MM:SS.FFF'))
rec_signal = double(rxLogAllData.Buffer)/2^15;
save('tx_backup.mat', 'rec_signal');
spectrogram(rec_signal,256,0,256,Fs,'yaxis','centered');
message_out = debug_LoRa_Rx(rec_signal,BW,SF,2,Fs, Fc-fc) ;
%% Message Out
disp(['Message Received = ' char(message_out)])