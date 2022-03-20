function m = k_pdf(x, beta, alpha)
m = 2/(alpha*gamma(beta+1)).*(x/(2*alpha)).^(beta+1).*besselk(beta, x/alpha);
m(find(isnan(m))) = 0;