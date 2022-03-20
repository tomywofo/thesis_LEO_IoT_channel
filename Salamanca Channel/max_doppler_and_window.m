function [tau, time_20, f_d_max_b, f_d_max_g] = max_doppler_and_window(theta_v, theta_max, height, inclination, fc)
% Compute Maximum Doppler value, along with visibility window duration.
% Receives orbital parameters as input, as well as theta_v, defined as the
% minimum inclination angle to establish a link with the satellite.
G = 6.67430e-11;
M_e = 5.9722e24;
R_e = 6371e3;
r = R_e+height;
omega_e = 2*pi/(24*3600);
omega_s = sqrt(G*M_e/r^3);
omega_f = omega_s-omega_e*cos(inclination*pi/180);
tau = 2/omega_f*acos((cos(acos(R_e/r*cos(theta_v*pi/180))-theta_v*pi/180))./(cos(acos(R_e/r*cos(theta_max*pi/180))-theta_max*pi/180)));
[~, f_d_max_b] = doppler_est(-tau/2, theta_max, height, inclination, fc);
time_20 = find_20(theta_max, height, inclination, tau);
[~, f_d_max_g] = doppler_est(time_20, theta_max, height, inclination, fc);