function k = rayleigh_lognormal_cdf_amplitude(x, mu, sigma)
scalar_fun = @(x) integral(@(y) rayleigh_lognormal_pdf_amplitude(y, mu, sigma), 0, x, 'RelTol', 1e-6,'ArrayValued', true);
cdf_func = @(x) arrayfun(scalar_fun, x);
k = cdf_func(x);