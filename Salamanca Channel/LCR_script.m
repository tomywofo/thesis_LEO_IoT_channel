SNR_db = -9:0.5:8;
SNR = 10.^(SNR_db/10);
theta_v = 10;
fc = 2.4e9;
theta_max = 90;
height = 400e3;
inclination = 0;

tau = max_doppler_and_window(theta_v, theta_max, height, inclination);
f_d_max = doppler_est(-tau/2, theta_max, height, inclination, fc);
Nn_k = LCR(f_d_max, SNR, 2, 10^(2/10));

semilogy(SNR_db, Nn_k/abs(max(Nn_k)))
axis([-40, 10, 10^-2, 10^0])