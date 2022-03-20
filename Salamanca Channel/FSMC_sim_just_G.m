clear all
k = 4;
fc = 2.4e9;
theta_v = 5;
theta_max = 90;
pi_k = 1/k;
Rs = 10^6; 
height = 680e3;
N_bits = 10^5;
montecarlo_N = 200;
inclination = 0;
min_SNR = -10;
step_SNR = 5;
max_SNR = 30;
m = [1, 2, 3, 4, 5, 10, 15];
SNR_vec = min_SNR:step_SNR:max_SNR;
BER_vec = zeros(length(SNR_vec),length(m));
[tau, time_20, f_d_b_max, f_d_g_max] = max_doppler_and_window(theta_v, theta_max, height, inclination, fc);

for i = 1:length(SNR_vec)
    SNR = 10^(SNR_vec(i)/10);
    for j = 1:length(m)   
        %thresholds = pdf_equipartition('Nakagami', k, [m(j), SNR]);
        thresholds = threshold_calc_g(k, m(j), SNR, 1e-10);     
        lcr_vals = LCR(f_d_g_max, thresholds, m(j), SNR);
        [previous, stay, next] = transition_probs(lcr_vals, Rs, pi_k);
        error_probs = error_prob_g(thresholds, m(j), SNR, sqrt(2), pi_k);
        total_BER = 0;
        for l = 1:montecarlo_N
            if mod(l, 50) == 0   
                fprintf(1, 'SNR = %d, m = %d, Montecarlo Nº: %d\n', SNR_vec(i), m(j), l)
            end
            total_BER = total_BER + markov_simulate(previous, stay, next, error_probs, N_bits);
        end
        BER_vec(i, j) = total_BER/montecarlo_N;
    end
end
figure()
for i = 1:length(m)    
    semilogy(SNR_vec, BER_vec(:, i), '*--', 'DisplayName', sprintf('Simulated Nakagami FSMC m = %d', m(i)))
    hold on
end
SNR_analytical = min_SNR:0.1:max_SNR;
SNR_db = 10.^(SNR_analytical/10);
analytical = 1/2*(1-sqrt(SNR_db./(1+SNR_db)));
plot(SNR_analytical, analytical, 'r-', 'DisplayName', 'Analytical Rayleigh')
ber_awgn = berawgn(SNR_analytical, 'psk', 2, 'nondiff');
plot(SNR_analytical, ber_awgn, 'k-', 'DisplayName', 'AWGN')
% SNR_analytical = -10:0.1:30;
% SNR_db = 10.^(SNR_analytical/10);
% analytical = 1/2*(1-sqrt(SNR_db./(1+SNR_db)));
% plot(SNR_analytical, analytical, 'r-', 'DisplayName', 'Analytical Rayleigh')
xlabel('SNR (dB)')
ylabel('BER')
legend()
% xlim([-10, 12])
ylim([10^(-4), 10^0])
title('Simulacion FSMC, Sector G')
