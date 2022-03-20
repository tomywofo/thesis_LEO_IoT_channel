function m = LCR_k_v2(f_d_max, thresholds, EsNo, beta_val)
power_tx = 1;
alpha_val = sqrt(power_tx)/2;
C = f_d_max*4*thresholds.^(beta_val/2)/((gamma(beta_val+1))^2)*(1/(4*alpha_val^2*EsNo))^(beta_val+2).*besselk(beta_val,2*(1/(4*alpha_val^2*EsNo))^(1/2)*thresholds.^(1/2));
D = 2*1/(4*alpha_val^2*EsNo)^(1/2);
m = 2*C/D^(beta_val+2)*2^(beta_val)*gamma(beta_val+1);
m(isnan(m)) = 0;