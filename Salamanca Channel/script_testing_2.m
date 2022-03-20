figure()
SNR_db = -3;
SNR = 10^(SNR_db/10);
%param_1 = 1.5;
thresholds = threshold_calc_b(k, SNR);
x = 0:0.05:5;
nakagami = @(z) pdf('Rayleigh', z, SNR);
y = nakagami(x);
plot(x,y, 'b', 'DisplayName', 'PDF')
title('Equipartición de distribución Rayleigh - PDF')
xlabel('x')
ylabel('PDF')
hold on
for i=1:length(thresholds)
    plot([thresholds(i,1),thresholds(i,1)], [0, nakagami(thresholds(i,1))], 'r-')
end
xlim([0.0, 3.0])

figure()
nakagami_cdf = @(z) cdf('Rayleigh', z, param_1);
y = nakagami_cdf(x);

plot(x,y, 'b', 'DisplayName', 'CDF')
xlabel('x')
ylabel('CDF')
title('Equipartición de distribución Rayleigh - CDF')
hold on
for i=1:length(thresholds)
    plot([thresholds(i,1),thresholds(i,1)], [0, nakagami_cdf(thresholds(i,1))], 'r-')
end