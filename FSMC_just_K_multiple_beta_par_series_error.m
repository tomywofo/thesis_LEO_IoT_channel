clear all
addpath('C:\Users\nunez\Desktop\Universidat\Magister\MATLAB\Salamanca Channel');
n_states = 4;
fc = 2.4e9;
theta_v = 5;
theta_max = 90;
pi_k = 1/n_states;
Rs = 10^6; 
height = 680e3;
N_bits = 10^6;
montecarlo_N = 400;
inclination = 0;
min_SNR = -10;
step_SNR = 5;
max_SNR = 30;
SNR_vec = min_SNR:step_SNR:max_SNR;
BER_vec = zeros(length(SNR_vec),1);
[tau, time_20, f_d_max_b, f_d_max_g] = max_doppler_and_window(theta_v, theta_max, height, inclination, fc);
% rec_power = 1;
% alpha = sqrt(rec_power)/2;
% lambda^2 = derivative(psi(1+beta)) = psi(1+beta, k)
beta_vals = [-0.37, 0.35, 2.001, 5.001];
style_sim = ["ko", "k*", "k+", "kx"];
style_an = ["r-", "r--", "r-.", "r:"];
style_gen = ["k-", "k--", "k-.", "k:"];
SNR_lin = 10.^(SNR_vec/10);
figure()
analytical_ray = 1/2*(1-sqrt(SNR_lin./(1+SNR_lin)));
semilogy(SNR_vec, analytical_ray', 'b-', 'DisplayName', 'Analytical Rayleigh')
hold on
for m = 1:length(beta_vals)
    beta_val = beta_vals(m);
    k_val = beta_val + 1;
    % beta = 1;
    parfor i = 1:length(SNR_vec)
        EsNo = 10^(SNR_vec(i)/10)/k_val;
        thresholds = threshold_calc_gen(n_states, 1e-12, @(z) k_SNR_cdf_an(z, EsNo, beta_val), max(3*EsNo*abs(beta_val),1));
        lcr_vals = LCR_k_v2(f_d_max_b, thresholds, EsNo, beta_val);
        [previous, stay, next] = transition_probs(lcr_vals, Rs, pi_k);
        error_probs = error_prob_K_series(thresholds, pi_k, EsNo, beta_val);
        total_BER = 0;
        for j = 1:montecarlo_N
            total_BER = total_BER + markov_simulate(previous, stay, next, error_probs, N_bits);
            if mod(j, 50) == 0   
                fprintf(1, 'SNR = %d, beta = %.2f, Montecarlo Nº: %d\n', SNR_vec(i), beta_val, j)
            end
        end
        BER_vec(i) = total_BER/montecarlo_N;
    end

    semilogy(SNR_vec, BER_vec', style_sim(m), 'DisplayName', ['Simulated FSMC K-fading (beta=', num2str(beta_val), ')'])
    analytical_gen = zeros(1, length(SNR_lin));
    for i = 1:length(analytical_gen)
        analytical_gen(i) = error_prob_gen_K_bpsk(SNR_lin(i)/k_val, beta_val);
    end
    analytical = error_prob_K_bpsk_v2(SNR_lin/k_val, beta_val);
    %semilogy(SNR_vec, analytical_gen', style_gen(m), 'DisplayName', ['Analytical K-fading Gen (beta=', num2str(beta_val), ')'])
    semilogy(SNR_vec, analytical', style_an(m), 'DisplayName', ['Analytical K-fading (beta=', num2str(beta_val), ')'])
    % analytical_ray = 1/2*(1-sqrt(((beta_val+1)*SNR_lin)./(1+((beta_val+1)*SNR_lin))));
    %SNR_lin = SNR_lin*(beta_val+1);

    % ber_ray = berfading(SNR_vec, 'psk', 2, 1);
    % plot(SNR_vec, ber_ray, 'k-')
    %beta = SNR/(4*alpha^2) - 1;
end
xlabel('SNR (dB)')
ylabel('BER')
legend()
title('FSMC Simulation for B-sector with K-distributed fading')
set(gcf, 'Position', [730 568 845 699])
grid on
save_name = 'ber_FSMC_k_fading_many_betas_par_series_error';
saveas(gcf, ['figuras\', save_name],'png')
saveas(gcf, ['figuras\', save_name],'epsc')
saveas(gcf, ['figuras\', save_name],'fig')