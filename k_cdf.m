function m = k_cdf(x, beta, alpha)
m = 1 - 2^(-beta)/gamma(beta+1)*(x/alpha).^(beta+1).*besselk(beta+1, x/alpha);
m(find(isnan(m))) = 0;