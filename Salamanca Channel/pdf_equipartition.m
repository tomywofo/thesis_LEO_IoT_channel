function thresholds = pdf_equipartition(pdf, k , params)
% PDF Equipartition:  Calculate thresholds by forcing equal steady state 
% probability. Only Rayleigh and Nakagami PDFs are accepted
%
% thresh = pdf_equipartition('Nakagami', 10, [5,3]) stores in thresh the 
% threshold values to be used for a 10-state FSMC, under a Nakagami PDF,
% with parameters mu = 5 and omega = 3.
%
% thresh = pdf_equipartition('Rayleigh', 8, 0.5) stores in thresh the 
% threshold values to be used for an 8-state FSMC, under a Rayleigh PDF,
% with parameter b = 0.5.

if strcmp(pdf, 'Nakagami')
    calc_cdf = @(x) cdf(pdf, x, params(1), params(2));
    calc_icdf = @(x) icdf(pdf, x, params(1), params(2));
elseif strcmp(pdf, 'Rayleigh')
    calc_cdf = @(x) cdf(pdf, x, params);
    calc_icdf = @(x) icdf(pdf, x, params);
else
    disp('Invalid pdf')
    return
end
thresholds = zeros(k+1,2);
current_val = 0;
current_cdf = 0;
target_cdf = 1/k;
for i=2:k
    next_val = calc_icdf(current_cdf+target_cdf);
    thresholds(i, 1) = current_val;
    thresholds(i, 2) = next_val;
    current_val = next_val;
    current_cdf = calc_cdf(current_val);
end
thresholds = thresholds(:,2);
thresholds(k+1)=Inf;