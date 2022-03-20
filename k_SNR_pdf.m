function m = k_SNR_pdf(gamma_val, EsNo, beta_val)
power_tx = 1;
alpha = sqrt(power_tx)/2;
% m = 2/gamma(beta+1)*(sqrt(gamma_val/EsNo)).^(beta+2)*(1/(2*alpha)).^(beta + 2).*besselk(beta, 1/alpha*sqrt(gamma_val/EsNo));
%m = 2*(gamma_val.^(beta/2))./((2*alpha)^(beta+2)*gamma(beta+1)*(sqrt(EsNo))^(beta+2)).*besselk(beta, 1/alpha*sqrt(gamma_val/EsNo));
%m = 2*(1/(EsNo*4*alpha^2))^((beta_val+2)/2)*gamma_val.^(beta_val/2)/(gamma(beta_val+1)).*besselk(beta_val, 2*(1/(EsNo*4*alpha^2))^(1/2)*gamma_val.^(1/2));
m = 2^(-beta_val-1)/(alpha^2*EsNo*gamma(beta_val+1))*(gamma_val/(alpha^2*EsNo)).^(beta_val/2).*besselk(beta_val, (gamma_val/(alpha^2*EsNo)).^(1/2));
m(gamma_val==0) = limit_at_0_k_pdf(EsNo, beta_val);