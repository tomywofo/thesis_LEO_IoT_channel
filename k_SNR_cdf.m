function k = k_SNR_cdf(x, EsNo, beta_val)
scalar_fun = @(x) integral(@(y) k_SNR_pdf(y, EsNo, beta_val), 0, x, 'RelTol', 1e-8,'ArrayValued', true);
cdf_func = @(x) arrayfun(scalar_fun, x);
k = cdf_func(x);
k(isnan(k)) = 0;