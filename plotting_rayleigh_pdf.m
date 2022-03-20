mu_db = 5;
sigma_db = 0.25;
mu_lin = 10^(mu_db/10);
sigma_lin = 10^(sigma_db/10);
SNR_low_db = -20;
SNR_high_db = 20;
SNR_db_linspaced = SNR_low_db:0.5:SNR_high_db;
SNR_lin = 10.^(SNR_db_linspaced);
SNR_lin_for_plot = -5:0.05:10;
pdf_lin = rayleigh_lognormal_pdf(SNR_lin_for_plot, mu_lin, sigma_lin);
plot(SNR_lin_for_plot, pdf_lin)
% mu = 1;
% SNR_db_linspaced
% SNR_lin_logspaced = logspace(log10(SNR_low_lin), log10(SNR_high_lin), 1000)
