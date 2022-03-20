function e_gk = error_prob_b(thresholds, SNR, a_mod, pi_k)
% Crossover Error Probability: Compute bit error probability for a FSMC
% Channel Model, based on the threshold values, SNR, modulation scheme and
% steady state probabilities.

K = length(thresholds)-1;
e_gk = zeros(K,1);
beta_k = exp(-thresholds/SNR).*(1-normcdf(a_mod*sqrt(thresholds))) + (sqrt(a_mod^2*SNR/(a_mod^2*SNR+2)))*normcdf(sqrt(thresholds*(a_mod^2*SNR+2)/SNR));
for i=1:K
    e_gk(i) = (beta_k(i) - beta_k(i+1))/pi_k;
end
