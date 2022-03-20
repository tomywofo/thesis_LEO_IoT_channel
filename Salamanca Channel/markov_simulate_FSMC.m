function [BER] = markov_simulate_FSMC(prev_g, stay_g, next_g, error_probs_g, prev_b, stay_b, next_b, error_probs_b, N_bits, good_bits_start, good_bits_end)
current_state = 1;
BER = 0;
prev = prev_b;
stay = stay_b;
next = next_b;
error_probs = error_probs_b;
for i=1:N_bits
    if i == good_bits_start
        prev = prev_g;
        stay = stay_g;
        next = next_g;
        error_probs = error_probs_g;
    elseif i == good_bits_end
        prev = prev_b;
        stay = stay_b;
        next = next_b;
        error_probs = error_probs_b;
    end
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