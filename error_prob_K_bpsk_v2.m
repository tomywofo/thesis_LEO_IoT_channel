function m = error_prob_K_bpsk_v2(EsNo, beta_val)
power_tx = 1;
alpha = sqrt(power_tx)/2;
m = 1/2*(1-(gamma(beta_val+3/2))./((4*alpha^2*EsNo).^(beta_val+1)* gamma(beta_val+1)).*kummerU(beta_val+3/2, beta_val+2, 1./(4*alpha^2*EsNo)));