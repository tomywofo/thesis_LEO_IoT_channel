function Nnk = LCR_k(f_d_max, thresholds, beta, alpha)
Nnk = 2*f_d_max*sqrt(pi)/(gamma(beta+1))^2.*(thresholds/(2*alpha)).^(beta+1).*besselk(beta, thresholds/alpha)*gamma((3+2*beta)/2);
Nnk(end) = 0;