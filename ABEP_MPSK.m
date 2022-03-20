function k = ABEP_MPSK(gamma_s, M, mu, sigma)
theta = 2/(log(4)/log(2));
lambda = 2*(sin(pi/M))^2;
k = integral(@(beta) theta*qfunc(sqrt(lambda*beta*gamma_s))*rayleigh_lognormal_pdf(beta, mu, sigma), 0, Inf, 'RelTol', 1e-6, 'ArrayValued', true);
