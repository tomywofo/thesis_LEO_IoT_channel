function e_gk = error_prob_g(thresholds, m, SNR, a_mod, pi_k)
% Crossover Error Probability: Compute bit error probability for a FSMC
% Channel Model, based on the threshold values, SNR, modulation scheme and
% steady state probabilities.

K = length(thresholds)-1;
e_gk = zeros(K,1);
Ik = zeros(K+1,1);
Fgamma = ((0:K)/K).';
for k = 1:K  
    Ii = 0:100;
    Ii = (-1).^(Ii).*((m/SNR)*(a_mod^2/2)).^(m+Ii)./(factorial(Ii).*(m+Ii)).*gammainc(a_mod^2/2*thresholds(k), m+Ii+1/2).*gamma(m+Ii+1/2);
    Ik(k) = 1/(2*sqrt(pi)*gamma(m))*sum(Ii);
end
beta_k = Fgamma.*(1-normcdf(a_mod*sqrt(thresholds)))+Ik;
for i=1:K-1
    e_gk(i) = (beta_k(i+1) - beta_k(i))/pi_k;
end
constant_C = (SNR*a_mod^2)/(2*m);
constant_value = sqrt(constant_C/(1 + constant_C));
sum_all_m = 0;
for k = 0:m-1
    sum_all_m = sum_all_m + nchoosek(2*k, k)*((1-constant_value^2)/4)^k;
end
average_error_prob = (1/2)*(1-constant_value^1*sum_all_m);
e_gk(K) = (average_error_prob - e_gk(1:K-1)'*repmat(pi_k, [K-1,1]))/pi_k;
% e_gk(e_gk > 1) = 1;
