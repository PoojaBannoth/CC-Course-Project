function [r_k_next] = rk_next(r_k)
% update the intermediate variable r
r_k_next = 1 + sqrt(1 + 4 * (r_k^2)) / 2;
end

