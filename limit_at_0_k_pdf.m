function k = limit_at_0_k_pdf(EsNo, beta_val)
power_tx = 1;
alpha = sqrt(power_tx)/2;
const_1 = 2^(-beta_val-1)/(gamma(beta_val+1)*alpha^2*EsNo);
if beta_val > 0
    if rem(beta_val,1) ~= 0
        limit_val = 2^(beta_val-1)*pi*csc(pi*beta_val)/(gamma(1-beta_val));
    else
        limit_val = 2^(beta_val-1)*factorial(beta_val-1);
    end
else
    limit_val = Inf;
end
k = limit_val*const_1;