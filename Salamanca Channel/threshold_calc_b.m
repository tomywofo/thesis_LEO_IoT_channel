function thresholds = threshold_calc_b(k, SNR)
thresholds = zeros(k+1,1);
for i=1:k-1
    thresholds(i+1) = -SNR*log(exp(-thresholds(i)/SNR)-1/k);
end
thresholds(k+1) = Inf;