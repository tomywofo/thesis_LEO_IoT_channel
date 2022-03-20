function k = nakagami_cdf(gamma_val, m, SNR)
k = gammainc(m/SNR*gamma_val, m);