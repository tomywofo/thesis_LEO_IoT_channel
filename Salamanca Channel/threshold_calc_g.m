function thresholds = threshold_calc_g(k, m, SNR, eps)
thresholds = zeros(k+1,1);
cdf_gamma = @(gamma_val) gammainc(m/SNR*gamma_val, m);
% for i=1:k-1
%     upper = thresholds(i) + max(10, 3*SNR*m);
%     lower = thresholds(i);
%     test_val = (upper+lower)/2;
%     estim = cdf_gamma(test_val) - cdf_gamma(thresholds(i));
%     dif = estim-1/k;
%     count = 0;
%     while abs(dif) > eps
%         if dif > 0
%             upper = (upper+lower)/2;
%         else
%             lower = (upper+lower)/2;
%         end
%         test_val = (upper+lower)/2;
%         estim = cdf_gamma(test_val) - cdf_gamma(thresholds(i));
%         dif = estim-1/k;
%         count = count + 1;
%     end
%     thresholds(i+1) = test_val;
% end
target_diff = 1/k;
for i=1:k-1
    upper = thresholds(i) + 3*SNR*m;
    lower = thresholds(i);
    while cdf_gamma(upper) - cdf_gamma(lower) < target_diff
        upper = upper + max(3*SNR*m, lower);
        disp(upper)
    end
    test_val = (upper+lower)/2;
    estim = cdf_gamma(test_val) - cdf_gamma(thresholds(i));
    dif = estim-target_diff;
    count = 0;
    while abs(dif) > eps
        if dif > 0
            upper = test_val;
        else
            lower = test_val;
        end
        test_val = (upper+lower)/2;
        estim = cdf_gamma(test_val) - cdf_gamma(thresholds(i));
        dif = estim-target_diff;
        count = count + 1;
    end
    thresholds(i+1) = test_val;
end
thresholds(k+1) = Inf;