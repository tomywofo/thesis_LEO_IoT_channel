clear
mu_db = 1;
sigma_db = 0.25;
mu = 10^(mu_db/10);
sigma = 10^(sigma_db/10);
fc = 2.4e9;
thresholds = threshold_calc_rl(10, mu, sigma, 1e-6);
x = 0:0.1:17;
nakagami = @(z) rayleigh_lognormal_pdf(z, mu, sigma);
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
nakagami_cdf_aux = @(z) rayleigh_lognormal_cdf(z, mu, sigma);
y = nakagami_cdf_aux(x);

plot(x,y, 'b', 'DisplayName', 'CDF')
xlabel('x')
ylabel('CDF')
title('Equipartición de distribución Nakagami - CDF')
hold on
for i=1:length(thresholds)
    plot([thresholds(i,1),thresholds(i,1)], [0, nakagami_cdf_aux(thresholds(i,1))], 'r-')
end

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

