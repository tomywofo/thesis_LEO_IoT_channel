function thresholds = threshold_calc_gen(k, eps, cdf_func, test_delta)
% FIX: Params for initial estimate
thresholds = zeros(k+1,1);
target_diff = 1/k;
for i=1:k-1
    upper = thresholds(i) + test_delta;
    lower = thresholds(i);
    while cdf_func(upper) - cdf_func(lower) < target_diff
        upper = upper + max(test_delta, lower);
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