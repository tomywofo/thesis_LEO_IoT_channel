clear
bad_SNR_db = -3;
k = 8;
fc = 2.4e9;
bad_SNR = 10^(bad_SNR_db/10);
theta_v = 10;
theta_max = 60;
pi_k = 1/k;
Rs = 10^5; 
height = 1000e3;
inclination = 0;
thresholds = threshold_calc_b(k, bad_SNR);
thresholds_2 = threshold_calc_g(k, 1, bad_SNR, 1e-8);

[tau, time_20, f_d_max_b, f_d_max_g] = max_doppler_and_window(theta_v, theta_max, height, inclination, fc);
lcr_vals = LCR(100, thresholds, 1, bad_SNR);
lcr_vals_2 = LCR(100, thresholds_2, 1, bad_SNR);
[previous, stay, next] = transition_probs(lcr_vals, Rs, pi_k);
[previous_2, stay_2, next_2] = transition_probs(lcr_vals_2, Rs, pi_k);

