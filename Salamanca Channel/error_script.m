clear
bad_SNR_db = 5;
k = 10;
fc = 2.4e9;
SNR = 10^(bad_SNR_db/10);
theta_v = 10;
theta_max = 60;
pi_k = 1/k;
Rs = 10e5; 
height = 1000e3;
a_mod = sqrt(2);
inclination = 45;
thresholds = pdf_equipartition('Rayleigh', k, SNR);

test_hola = @(thresholds) exp(-thresholds/SNR).*(1-normcdf(a_mod*sqrt(thresholds))) + (sqrt(a_mod^2*SNR/(a_mod^2*SNR+2)))*normcdf(sqrt(thresholds*(a_mod^2+SNR+2)/SNR));
a = 0:0.05:10;
SNR_lala = test_hola(a);
plot(a, SNR_lala)