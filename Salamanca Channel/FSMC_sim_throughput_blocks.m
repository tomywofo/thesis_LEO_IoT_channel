clear all
k = 4;
fc = 2.4e9;
theta_v = 0;
theta_max = 90;
pi_k = 1/k;
Rs = 250e3; 
height = 680e3;
% N_bits = 10^4;
block_size = 64;
montecarlo_N = 200;
inclination = 0;
min_SNR = -2;
step_SNR = 1;
max_SNR = 12;
m = [1, 2, 3, 4, 5, 10, 15];
SNR_vec = min_SNR:step_SNR:max_SNR;
BER_vec = zeros(length(SNR_vec),length(m));
[tau, time_20, f_d_b_max, f_d_g_max] = max_doppler_and_window(theta_v, theta_max, height, inclination, fc);
N_blocks = round(tau/2*Rs/64);
good_block_start = round((tau/2-abs(time_20))*Rs/block_size);
just_b = zeros(length(SNR_vec), 1);
just_g = zeros(length(SNR_vec), 1);
for i = 1:length(SNR_vec)
    SNR = 10^(SNR_vec(i)/10);
    for j = 1:length(m)
        initial_time = -tau/2;
        %thresholds = pdf_equipartition('Nakagami', k, [m(j), SNR]);
        thresholds_g = threshold_calc_g(k, m(j), SNR, 1e-8);
        lcr_vals_g = LCR(f_d_g_max, thresholds_g, m(j), SNR);
        [previous_g, stay_g, next_g] = transition_probs(lcr_vals_g, Rs, pi_k);
        error_probs_g = error_prob_g(thresholds_g, m(j), SNR, sqrt(2), pi_k);
        thresholds_b = threshold_calc_b(k, SNR);
        lcr_vals_b = LCR(f_d_b_max, thresholds_b, 1, SNR);
        [previous_b, stay_b, next_b] = transition_probs(lcr_vals_b, Rs, pi_k);
        error_probs_b = error_prob_b(thresholds_b, SNR, sqrt(2), pi_k);
        total_BER = 0;
        for l = 1:montecarlo_N
            if mod(l, 50) == 0   
                fprintf(1, 'SNR = %d, m = %d, Montecarlo Nº: %d\n', SNR_vec(i), m(j), l)
            end
            total_BER = total_BER + markov_simulate_FSMC_blocks(previous_g, stay_g, next_g, error_probs_g, previous_b, stay_b, next_b, error_probs_b, N_blocks, block_size, good_block_start);
        end
        BER_vec(i, j) = total_BER/montecarlo_N;
        if m(j) == 1
            bad_BER = 0;
            for l = 1:montecarlo_N
                bad_BER = bad_BER + markov_simulate_blocks(previous_b, stay_b, next_b, error_probs_b, N_blocks, block_size);
            end
            just_b(i) = bad_BER/montecarlo_N;
        elseif m(j) == 15
            good_BER = 0;
            for l = 1:montecarlo_N
                good_BER = good_BER + markov_simulate_blocks(previous_g, stay_g, next_g, error_probs_g, N_blocks, block_size);
            end
            just_g(i) = good_BER/montecarlo_N;
        end
    end
end
analytical_SNR = min_SNR:0.1:max_SNR;
awgn_channel = berawgn(analytical_SNR, 'psk', 2, 'nondiff');
figure()
for i = 1:length(m)    
    semilogy(SNR_vec, BER_vec(:, i), '*-', 'DisplayName', sprintf('Simulated Nakagami FSMC m = %d', m(i)))
    hold on
end
% SNR_analytical = -10:0.1:30;
% SNR_db = 10.^(SNR_analytical/10);
% analytical = 1/2*(1-sqrt(SNR_db./(1+SNR_db)));
% plot(SNR_analytical, analytical, 'r-', 'DisplayName', 'Analytical Rayleigh')
xlabel('SNR (dB)')
ylabel('BER')
legend()
title('Simulacion FSMC - BER')

figure()
hold on
for i = 2:length(m)
    throughput_norm = throughput(block_size, block_size, BER_vec(:,i));
    plot(SNR_vec, throughput_norm, '*-', 'DisplayName', sprintf('Simulated Nakagami FSMC m = %d', m(i)))
end
bad_throughput = throughput(block_size, block_size, just_b);
good_throughput = throughput(block_size, block_size, just_g);
awgn_throughput = throughput(block_size, block_size, awgn_channel);
plot(SNR_vec, bad_throughput, '*-', 'DisplayName', sprintf('Only B sector simulation'))
plot(SNR_vec, good_throughput, '*-', 'DisplayName', sprintf('Only G sector simulation'))
plot(analytical_SNR, awgn_throughput, '-', 'DisplayName', sprintf('AWGN analytical'))
xlabel('SNR (dB)')
ylabel('Throughput (eta)')
legend()
title('Throughput teórico en FSMC')    