function k = rayleigh_lognormal_pdf_amplitude(x, mu, sigma)
k = integral(@(gamma) x'./gamma.*exp(-(x.^2)'./(2*gamma)).*1./(sqrt(2*pi)*sigma*gamma).*exp(-(log(gamma)-mu).^2./(2*sigma.^2)), 0, Inf, 'RelTol', 1e-6, 'ArrayValued', true);

