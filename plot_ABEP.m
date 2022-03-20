M = 4; % QPSK
mu_db = 1;
sigma_db = 0.5;
mu_lin = 10^(mu_db/10);
sigma_lin = 10^(sigma_db/10);
SNR_high_db = 40;
SNR_db_linspaced = 0:0.5:SNR_high_db;
SNR_lin = 10.^(SNR_db_linspaced/10);
ABEP = ABEP_MPSK(SNR_lin, M, mu_lin, sigma_lin);
semilogy(SNR_db_linspaced, ABEP, '*')
xlim([0, 40]);
ylim([10^-5, 1])