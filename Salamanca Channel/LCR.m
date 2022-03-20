function Nnk = LCR(f_d_max, thresholds, m, SNR)
Nnk = sqrt(2*pi)/gamma(m)*f_d_max*(m/SNR*thresholds).^(m-1/2).*exp(-m/SNR*thresholds);
Nnk(end) = 0;