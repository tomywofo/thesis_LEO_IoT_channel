clear all
addpath('C:\Users\nunez\Desktop\Universidat\Magister\MATLAB\Salamanca Channel');

beta_1 = -0.65;
beta_2 = -0.37;
beta_3 = 2.5;
k = 4;
pi_k = 1/k;
fc = 2.4e9;
EsNo_db = 30;
EsNo = 10^(EsNo_db/10);
thresholds = threshold_calc_gen(k, 1e-12, @(z) k_SNR_cdf_an(z, EsNo, beta_3), max(3*EsNo*abs(beta_3),1));
%x = 0:1e-4:3;
x = linspace(0,ceil(thresholds(end-1)), 1e4);
% SNR_db = 0:0.05:30;
% SNR = 10.^(SNR_db/10);
f_d_max = 80;
Rs = 10^5;
k_dist = @(z) k_SNR_pdf(z, EsNo, beta_3);
y = k_dist(x);
figure(1)
plot(x,y, 'b', 'DisplayName', 'pdf')
title('Pdf curve of a K-distributed SNR')
xlabel('x')
ylabel('pdf')
xlim([0, ceil(thresholds(end-1))])
ylim([0, ceil(max(k_dist(thresholds(2:end-1))))])
hold on
for i=1:length(thresholds)
    plot([thresholds(i),thresholds(i)], [0, k_dist(thresholds(i))], 'r-')
end
set(gcf, 'Position', [730 568 845 699])
grid on
saveas(gcf, 'k_pdf_eq.eps', 'epsc')

figure(2)
k_distr_aux =  @(z) k_SNR_cdf_an(z, EsNo, beta_3);
y = k_distr_aux(x);

plot(x,y, 'b', 'DisplayName', 'cdf')
xlabel('x')
ylabel('cdf')
title('Equipartition of cdf for K-distributed SNR')
hold on
for i=1:length(thresholds)
    plot([thresholds(i,1),thresholds(i,1)], [0, k_distr_aux(thresholds(i,1))], 'r-')
end
set(gcf, 'Position', [730 568 845 699])
grid on
xlim([0, ceil(thresholds(end-1))])
ylim([0, 1])
saveas(gcf, 'k_cdf_eq.eps', 'epsc')

figure(3)
k_lcr = LCR_k_v2(f_d_max, thresholds, EsNo, beta_3);
plot(1:k, k_lcr(2:end))
[previous, stay, next] = transition_probs(k_lcr, Rs, 1/k);

bit_error_rate_1 = error_prob_K_series(thresholds, pi_k, EsNo, beta_1);
bit_error_rate_2 = error_prob_K_series(thresholds, pi_k, EsNo, beta_2);
bit_error_rate_3 = error_prob_K_series(thresholds, pi_k, EsNo, beta_3);
bit_error_rate_3_v2 = error_prob_K(thresholds, pi_k, EsNo, beta_3);

% 
% plot(SNR_db, log10(bit_error_rate_1), 'b--','DisplayName', 'Beta = -0.65')
% hold on
% plot(SNR_db, log10(bit_error_rate_2), 'b-','DisplayName', 'Beta = -0.37')
% plot(SNR_db, log10(bit_error_rate_3), 'b-.','DisplayName', 'Beta = 0.35')
% ylim([-4, -0.5])
% ylabel('log(Average MSK BER)')
% xlabel('Eb/N0 [dB]')
% title('Average BER for a K-fading channel using MSK')
% legend('show')



% figure()
% theta1 = 20;
% theta2 = 90;
% time = -8*60:1:8*60;
% doppler_estim = doppler_est(time, 20, 400e3, 15, fc );
% doppler_estim_2 = doppler_est(time, 90, 400e3, 15, fc);
% plot(time/60, doppler_estim, 'DisplayName', sprintf('Theta max = %dº', theta1))
% hold on
% plot(time/60, doppler_estim_2, 'r-', 'DisplayName', sprintf('Theta max = %dº', theta2))
% title('Estimación efecto Doppler - Órbita de 400 Km')
% legend()
% xlabel('Tiempo (min)')
% ylabel('Delta f/f')

% figure()
% theta_max = 10:1:90;
% visibility = max_doppler_and_window(10, theta_max, 1000e3, 15, fc);
% plot(theta_max, visibility/60)
% title('Ventana de visibilidad en función de theta max (órbita de 1000 km)')
% xlabel('Theta max (deg)')
% ylabel('Ventana de visibilidad (min)')
