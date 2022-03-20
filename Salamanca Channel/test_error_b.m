clear
bad_SNR_db = -3;
k = 4;
fc = 2.4e9;
bad_SNR = 10^(bad_SNR_db/10);
theta_v = 10;
theta_max = 60;
pi_k = 1/k;
Rs = 10e5; 
height = 1000e3;
inclination = 45;
thresholds = pdf_equipartition('Rayleigh', k, bad_SNR);

tau = max_doppler_and_window(theta_v, theta_max, height, inclination);
[f_delta_max, f_d_max] = doppler_est(-tau/2, theta_max, height, inclination, fc);

lcr_vals = LCR(f_d_max, thresholds, 1, bad_SNR);
[previous, stay, next] = transition_probs(lcr_vals, Rs, pi_k);

error_vals = error_prob_b(thresholds, bad_SNR, sqrt(2), pi_k);