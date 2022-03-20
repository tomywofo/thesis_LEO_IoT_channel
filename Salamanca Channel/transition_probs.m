function [previous, stay, next]= transition_probs(LCR_vals, Rs, pi_k)
% Compute transition probabilities through formulas from FSMC-Channel paper
% (Salamanca et al, 2019).

K = length(LCR_vals)-1;
previous = zeros(K-1, 1);
next = zeros(K-1, 1);
for i=1:K-1
    previous(i) = LCR_vals(i+1)./(Rs*pi_k);
    next(i) = LCR_vals(i+1)./(Rs*pi_k);
end
previous = [0; previous];
next = [next; 0];
stay = 1 - next - previous;