function [BER] = markov_simulate_blocks(prev, stay, next, error_probs, N_blocks, block_size)
current_state = 1;
BER = 0;
for i=1:N_blocks
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