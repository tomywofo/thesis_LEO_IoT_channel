clear
good_SNR_db = -10;
k = 4;
fc = 2.4e9;
m = 3;
good_SNR = 10^(good_SNR_db/10);
theta_v = 5;
theta_max = 90;
pi_k = 1/k;
Rs = 10^5; 
height = 1000e3;
inclination = 45;
thresholds = threshold_calc_g(k, m, good_SNR, 1e-8);

tau = max_doppler_and_window(theta_v, theta_max, height, inclination);
[f_delta_max, f_d_max] = doppler_est(-tau/2, theta_max, height, inclination, fc);

lcr_vals = LCR(f_d_max, thresholds, 1, good_SNR);
[previous, stay, next] = transition_probs(lcr_vals, Rs, pi_k);

error_vals = error_prob_g(thresholds, m, good_SNR, sqrt(2), pi_k);