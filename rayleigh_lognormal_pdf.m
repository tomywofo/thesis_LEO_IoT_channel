function k = rayleigh_lognormal_pdf(x, mu, sigma)
k = integral(@(gamma) 1./gamma.*exp(-x'./gamma).*1./(sqrt(2*pi)*sigma*gamma).*exp(-(log(gamma)-mu).^2./(2*sigma.^2)), 0, Inf, 'RelTol', 1e-6, 'ArrayValued', true);

