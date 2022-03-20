function [BER] = markov_simulate(prev, stay, next, error_probs, N_bits)
current_state = 1;
BER = 0;
for i=1:N_bits
    error_gen = rand();
    if error_gen <= error_probs(current_state)
        BER = BER + 1;
    end
    trans_gen = rand();
    if trans_gen >= next(current_state) + stay(current_state)
        current_state = current_state - 1;
    elseif trans_gen >= stay(current_state)
        current_state = current_state + 1;
    else
        current_state = current_state;
    end
end
BER = BER/N_bits;