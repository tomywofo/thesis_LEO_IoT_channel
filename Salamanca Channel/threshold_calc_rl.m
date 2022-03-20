function thresholds = threshold_calc_rl(k, mu, sigma, eps)
thresholds = zeros(k+1,1);
cdf_func = @(z) rayleigh_lognormal_cdf(z, mu, sigma);
target_diff = 1/k;
for i=1:k-1
    upper = thresholds(i) + 3*mu*sigma;
    lower = thresholds(i);
    while cdf_func(upper) - cdf_func(lower) < target_diff
        upper = upper + max(3*mu*sigma, lower);
    end
    test_val = (upper+lower)/2;
    estim = cdf_func(test_val) - cdf_func(thresholds(i));
    dif = estim-target_diff;
    count = 0;
    while abs(dif) > eps
        if dif > 0
            upper = test_val;
        else
            lower = test_val;
        end
        test_val = (upper+lower)/2;
        estim = cdf_func(test_val) - cdf_func(thresholds(i));
        dif = estim-target_diff;
        count = count + 1;
    end
    thresholds(i+1) = test_val;
end
thresholds(k+1) = Inf;