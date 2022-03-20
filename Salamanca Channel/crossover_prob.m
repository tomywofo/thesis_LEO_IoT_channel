%% Crossover probability

a_parameter_bpsk_coherent = sqrt(2);
naka.beta_k = zeros(naka.numberStates,1);
for k=1:naka.numberStates   
    naka.gamma_cdf(k) = ...
        (1/gamma(naka.m_factor))*...
        function_gamma_incomplet_integral(naka.m_factor,naka.m_factor*naka.gamma_ep(k)/naka.gamma_average);    
    gammagamma_incomplet = 0;
    for n=0:numInf
        incomplet_gamma_integral = function_gamma_incomplet_integral(naka.m_factor+n+1/2,...
            (a_parameter_bpsk_coherent^2/2)*naka.gamma_ep(k));
        gammagamma_incomplet = gammagamma_incomplet + ...
            (((-1)^n*((naka.m_factor/naka.gamma_average)*(a_parameter_bpsk_coherent^2/2))^(naka.m_factor+n))/...
            (factorial(n)*(naka.m_factor+n)))*incomplet_gamma_integral;
    end
    naka.I_k(k,1) = (1/(2*sqrt(pi)*gamma(naka.m_factor)))*gammagamma_incomplet;
    naka.beta_k(k) = naka.gamma_cdf(k)*qfunc(a_parameter_bpsk_coherent*sqrt(naka.gamma_ep(k))) + naka.I_k(k,1);
end

naka.crossoverProb = zeros(naka.numberStates,1);           
for k=1:naka.numberStates-1
    naka.crossoverProb(k) =(naka.beta_k(k+1)-naka.beta_k(k))/naka.steady_state_prob(k);
end

% Crossover probability for the last state
%
 
naka.constantC          = (naka.gamma_average*a_parameter_bpsk_coherent^2)/(2*naka.m_factor);
naka.constantValue      = sqrt(naka.constantC/(1+naka.constantC));
naka.average_error_prob = 0;
naka.sumAll_m           = 0;
for k=0:naka.m_factor-1
        naka.sumAll_m = naka.sumAll_m + nchoosek(2*k,k)*((1-naka.constantValue^2)/4)^k;
end
naka.average_error_prob = (1/2)*(1-naka.constantValue^1*naka.sumAll_m);

naka.crossoverProb(end) = (naka.average_error_prob -...
    naka.crossoverProb(1:naka.numberStates-1)'*naka.steady_state_prob(1:naka.numberStates-1))   /...
    naka.steady_state_prob(end);

