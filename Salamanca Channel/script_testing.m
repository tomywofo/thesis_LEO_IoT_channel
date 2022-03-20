clear
SNR_db = -3;
SNR = 10^(SNR_db/10);
param_1 = 3;
param_2 = SNR;
fc = 2.4e9;
thresholds = threshold_calc_g(10, param_1,param_2, 1e-8);
x = 0:0.05:1.5;
nakagami = @(z) pdf('Nakagami', z, param_1, param_2);
y = nakagami(x);
figure()
plot(x,y, 'b', 'DisplayName', 'PDF')
title('Equipartición de distribución Nakagami - PDF')
xlabel('x')
ylabel('PDF')
hold on
for i=1:length(thresholds)
    plot([thresholds(i,1),thresholds(i,1)], [0, nakagami(thresholds(i,1))], 'r-')
end

figure()
nakagami_cdf_aux = @(z) nakagami_cdf(z, param_1, param_2);
y = nakagami_cdf_aux(x);

plot(x,y, 'b', 'DisplayName', 'CDF')
xlabel('x')
ylabel('CDF')
title('Equipartición de distribución Nakagami - CDF')
hold on
for i=1:length(thresholds)
    plot([thresholds(i,1),thresholds(i,1)], [0, nakagami_cdf_aux(thresholds(i,1))], 'r-')
end

figure()
theta1 = 20;
theta2 = 90;
time = -8*60:1:8*60;
doppler_estim = doppler_est(time, 20, 400e3, 15, fc );
doppler_estim_2 = doppler_est(time, 90, 400e3, 15, fc);
plot(time/60, doppler_estim, 'DisplayName', sprintf('Theta max = %dº', theta1))
hold on
plot(time/60, doppler_estim_2, 'r-', 'DisplayName', sprintf('Theta max = %dº', theta2))
title('Estimación efecto Doppler - Órbita de 400 Km')
legend()
xlabel('Tiempo (min)')
ylabel('Delta f/f')

% figure()
% theta_max = 10:1:90;
% visibility = max_doppler_and_window(10, theta_max, 1000e3, 15, fc);
% plot(theta_max, visibility/60)
% title('Ventana de visibilidad en función de theta max (órbita de 1000 km)')
% xlabel('Theta max (deg)')
% ylabel('Ventana de visibilidad (min)')

