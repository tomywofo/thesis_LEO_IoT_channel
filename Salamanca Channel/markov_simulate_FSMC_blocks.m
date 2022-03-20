function [BER] = markov_simulate_FSMC_blocks(prev_g, stay_g, next_g, error_probs_g, prev_b, stay_b, next_b, error_probs_b, N_blocks, block_size, good_block_start)
current_state = 1;
BER = 0;
prev = prev_b;
stay = stay_b;
next = next_b;
error_probs = error_probs_b;
for i=1:N_blocks
    if i == good_block_start
        prev = prev_g;
        stay = stay_g;
        next = next_g;
        error_probs = error_probs_g;
    end
    error_gen = rand(block_size, 1);
    BER = BER + length(error_gen(error_gen < error_probs(current_state)));
    trans_gen = rand();
    if trans_gen >= next(current_state) + stay(current_state)
        current_state = current_state - 1;
    elseif trans_gen >= stay(current_state)
        current_state = current_state + 1;
    else
        current_state = current_state;
    end
end
BER = BER/(N_blocks*block_size);