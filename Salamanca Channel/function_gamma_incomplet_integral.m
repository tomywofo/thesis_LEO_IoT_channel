function [ incomplet_gamma_integral ] = function_gamma_incomplet_integral(alpha,x)
% Incomplete Gamma Function Integral
%
% Inputs
% alpha:    
% x:        Integral limit
% Outputs
% incomplet_gamma_integral: Integral result fot incomplete gamma function  
% 
%
%   Incomplete gamma function definition [Ryzhik1980, eq.8.350(1), pg. 899]:
%       incomplete_gamma(alpha,x) = integral(exp(-t)*t^(alpha-1),0,x)
%       with Re{alpha} > 0


incomplete_gamma = @(t) exp(-t).*(t).^(alpha-1);        

incomplet_gamma_integral = integral(incomplete_gamma,0,x);


end

