addpath C:\Users\nunez\Desktop\Universidat\Magister\MATLAB\LoRa
addpath C:\Users\nunez\Desktop\Universidat\Magister\MATLAB\Patzold

%% Transmission Parameters
SF = 7 ;
BW = 125e3 ;
fc = 917.8e6 ;
Power = 0 ;
SNR = 30;
message = '';
for i = 1:15
    message = [message, 'Hello World ',  num2str(i), '\n'];
end
%% Channel Parameters
kappa_0 = 0.5;
sigma_0 = 0.86;
s_param = 0.5;
m = -0.25;
psi_0 = 2/pi*sigma_0^2*asin(kappa_0);
f_max = 2e6;
f_p = -1e6;
kappa_c = 50;
N1 = 10;
N2 = 10;

%% Sampling and carrier parameters
Fs = 1e6 ;
Fc = 917.8e6 ;