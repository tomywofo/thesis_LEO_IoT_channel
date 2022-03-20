function k = error_prob_K(thresholds, pi_k, EsNo, beta_val)
thresholds_up = thresholds(2:end);
thresholds_low = thresholds(1:end-1);
k = zeros(1,length(thresholds));
for i = 1:length(thresholds_up)
    k(i) = 1./pi_k.*integral(@(gamma_val) k_SNR_pdf(gamma_val, EsNo, beta_val).*qfunc(sqrt(2*gamma_val)), thresholds_low(i), thresholds_up(i), 'RelTol', 1e-12);
end
