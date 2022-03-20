function [f_delta, f_d_hz] = doppler_est(time, theta_max, height, inclination, fc)
% Doppler: Computes frequency delta due to Doppler effect. Receives orbital
% parameters and time (can be vector) as input. Time t = 0 corresponds to
% maximum elevation angle.
c = 3e8;
G = 6.67430e-11;
M_e = 5.9722e24;
R_e = 6371e3;
r = R_e+height;
omega_e = 2*pi/(24*3600);
omega_s = sqrt(G*M_e/r^3);
omega_f = omega_s-omega_e*cos(inclination*pi/180);
aux = cos(acos(R_e/r*cos(theta_max*pi/180))-theta_max*pi/180);
f_delta = -1/c*(R_e*r*sin(omega_f*time)*aux'*omega_f)./sqrt(R_e^2+r^2-2*R_e*r*cos(omega_f*time)*aux');
f_d_hz = f_delta*fc;