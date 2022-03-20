function e_gk = error_prob_K_series(thresholds, pi_k, EsNo, beta_val)
% Crossover Error Probability: Compute bit error probability for a FSMC
% Channel Model, based on the threshold values, SNR, modulation scheme and
% steady state probabilities.

K = length(thresholds)-1;
e_gk = zeros(K,1);
Jk = zeros(K+1,1);
Ik = zeros(K+1,1);
Fgamma = ((0:K)/K).';
power_tx = 1;
alpha = sqrt(power_tx)/2;
power_ct = alpha^2*EsNo;
sum_count = 100;
for h = 1:length(thresholds)
    k = 0:(sum_count-1);
    series_1 = gammainc(thresholds(h), (2*k+1)/2).*gamma((2*k+1)/2)./(gamma(k-beta_val).*factorial(k).*2.^(2*k-beta_val-1).*power_ct.^((2*k-1)/2));
    series_2 = gammainc(thresholds(h), (2*k+2*beta_val+3)/2).*gamma((2*k+2*beta_val+3) /2)./(gamma(k+beta_val+2).*factorial(k).*2.^(2*k+beta_val+1).*power_ct.^((2*k+2*beta_val+1)/2));
    %series_1 = (-1).^(series_1).*((m/SNR)*(a_mod^2/2)).^(m+series_1)./(factorial(series_1).*(m+series_1)).*gammainc(a_mod^2/2*thresholds(h), m+series_1+1/2).*gamma(m+series_1+1/2);
    Jk(h) = 2^(-beta_val-1)*pi*csc(pi*(beta_val+1))/(sqrt(2*pi)*gamma(beta_val+1)*sqrt(2*alpha^2*EsNo))*(sum(series_1)-sum(series_2));
    Ik = 0.5-qfunc(sqrt(2*thresholds))-Jk;
end
beta_k = Fgamma.*qfunc(sqrt(2*thresholds))+Ik;
for i=1:K
    e_gk(i) = (beta_k(i+1) - beta_k(i))/pi_k;
end
% constant_C = (SNR*a_mod^2)/(2*m);
% constant_value = sqrt(constant_C/(1 + constant_C));
% sum_all_m = 0;
% for h = 0:m-1
%     sum_all_m = sum_all_m + nchoosek(2*h, h)*((1-constant_value^2)/4)^h;
% end

% k_val = beta_val + 1;
average_error_prob = error_prob_K_bpsk_v2(EsNo, beta_val);
% average_error_prob = error_prob_gen_K_bpsk(EsNo*k_val, beta_val);
e_gk(K) = (average_error_prob - e_gk(1:K-1)'*repmat(pi_k, [K-1,1]))/pi_k;
e_gk(e_gk > 0.5) = 0.5;
