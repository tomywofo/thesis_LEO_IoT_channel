clear all
k = 4;
fc = 2.4e9;
theta_v = 5;
theta_max = 90;
pi_k = 1/k;
Rs = 10^6; 
height = 680e3;
N_bits = 10^4;
montecarlo_N = 200;
inclination = 0;
min_SNR = -10;
step_SNR = 5;
max_SNR = 20;
SNR_vec = min_SNR:step_SNR:max_SNR;
BER_vec = zeros(length(SNR_vec),1);
[tau, time_20, f_d_max_b, f_d_max_g] = max_doppler_and_window(theta_v, theta_max, height, inclination, fc);

for i = 1:length(SNR_vec)
    SNR = 10^(SNR_vec(i)/10);
    thresholds = threshold_calc_b(k, SNR);
    lcr_vals = LCR(f_d_max_b, thresholds, 1, SNR);
    [previous, stay, next] = transition_probs(lcr_vals, Rs, pi_k);
    error_probs = error_prob_b(thresholds, SNR, sqrt(2), pi_k);
    total_BER = 0;
    for j = 1:montecarlo_N
        total_BER = total_BER + markov_simulate(previous, stay, next, error_probs, N_bits);
    end
    BER_vec(i) = total_BER/montecarlo_N;
end

semilogy(SNR_vec, BER_vec, 'k*', 'DisplayName', 'Simulated FSMC B-sector')
hold on
SNR_analytical = min_SNR:0.1:max_SNR;
SNR_db = 10.^(SNR_analytical/10);
analytical = 1/2*(1-sqrt(SNR_db./(1+SNR_db)));
plot(SNR_analytical, analytical, 'b-', 'DisplayName', 'Analytical Rayleigh')
xlabel('SNR (dB)')
ylabel('BER')
legend()
title('Simulacion FSMC, Sector B')
    