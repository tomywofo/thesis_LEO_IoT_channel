function k = k_SNR_cdf_an(gamma_val, EsNo, beta)
power_tx = 1;
alpha_val = sqrt(power_tx)/2;
alpha_bar = alpha_val;
k = 1-(2^(-beta)/(gamma(beta+1)))*(gamma_val/(EsNo*alpha_bar^2)).^((beta+1)/2).*besselk(beta+1,sqrt(gamma_val/EsNo)/alpha_bar);
k(gamma_val==0) = 0;