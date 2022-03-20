close all

gamma_val = 0:0.01:6;
EsNo = 1;
beta_val = -0.7;
SNR_pdf = k_SNR_pdf(gamma_val, EsNo, beta_val);
SNR_cdf = k_SNR_cdf(gamma_val, EsNo, beta_val);
SNR_cdf_an = k_SNR_cdf_an(gamma_val, EsNo, beta_val);
figure()
plot(gamma_val, SNR_pdf)
title('PDF of received SNR under K-distribution')
xlabel('SNR [dB]')
ylabel('PDF')

figure()
hold on
plot(gamma_val, SNR_cdf, 'b-', 'DisplayName', 'Numerical')
plot(gamma_val, SNR_cdf_an, 'r', 'DisplayName', 'Analytical')
legend('show')
title('CDF of received SNR under K-distribution')
xlabel('SNR [dB]')
ylabel('CDF')
figure()
plot(gamma_val, abs(SNR_cdf-SNR_cdf_an))