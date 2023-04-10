function [y_out] = y(N_o, B, x)

y_out = N_o * B * (2^(x / B) - 1);
end

