B = 50 * 10^6; % in Mbps, Cannel Bandwidth
c_i = 1000; % cycles per bit
w_i = 500 * 10^3; % task in bits
f_c = 2 * 10^8; % cycles per second
f_i_l = 2 * 10^6; % local computing ability in cycles per second
N_o = 10^(-10); % in Watt
h_i_sq = 1;
E_i = 1.5; % in Joules, desired energy consumption
T_i = 600 * 10^(-3); % in seconds, desired latency
p_i_c = 10^(-3); % in Watt
m_i = 10^(-8); % in J/cycle
a_i = 0.5; % offloading ratio
max_iterations = 50;

R_i = B * log2(1 + p_i_c * h_i_sq / (N_o * B));
t_i = (w_i / R_i) * a_i;
d_i = w_i * a_i;


x = 1:max_iterations;
a = zeros(1,max_iterations);
R = zeros(1,max_iterations);
t = zeros(1,max_iterations);
p = zeros(1, max_iterations);
a(1, 1) = a_i;
R(1, 1) = R_i;
t(1, 1) = t_i;
disp(a_i);
disp(R_i);
disp(t_i);
disp(T_i);

t_i_c = w_i / R_i + c_i * w_i / f_c;
t_i_l = c_i * w_i / f_i_l;
energy_err = t_i_l + t_i_c;

lambda = zeros(1, max_iterations+1);
mu = zeros(1, max_iterations+1);

lambda(1, 1) = 100;
mu(1, 1) = 100;
lambda(1, 2) = 100;
mu(1, 2) = 100;
r = 2; % iteration index
rho = 5;
k = 0;
rk = zeros(1, max_iterations);
rk(1, 1) = 0;
rk(1, 2) = 0;

while r <= max_iterations
    % update r
    rk(1, r+1) = rk_next(rk(1, r));
    % update delta
    delta = lambda(1, r) + ((rk(1, r) - 1) / rk(1, r+1)) * (lambda(1, r) - lambda(1, r-1));
    % update epsilon
    epsilon = mu(1, r) + ((rk(1, r) - 1) / rk(1, r+1)) * (mu(1, r) - mu(1, r-1));

    % update lambda
    lambda(1, r+1) = delta + rho * (c_i * w_i * (1 - a_i) / f_i_l + ((w_i / R_i) + (c_i * w_i / f_c)) * a_i - T_i);
    
    %update mu
    mu(1, r+1) = ( epsilon + rho * (t_i / h_i_sq * y(N_o, B, d_i/t_i) + c_i * m_i * w_i * (1 - a_i) + energy_err - E_i));

    R(1, r) = B/log(2) * lambertw(0, (h_i_sq * (lambda(1, r) + (1+ mu(1, r)) * p_i_c)) / (N_o * B * exp(1) * (1+ mu(1, r))) - 1/exp(1)) + 1;
    a(1, r) = ((T_i * f_i_l) - (c_i * w_i)) / ((w_i * f_i_l / R_i) + (c_i * w_i * f_i_l / f_c) - (c_i * w_i));
    t(1, r) = w_i * a_i / R_i;
    p(1, r) = t(1, r) * y(N_o, B, d_i / t(1, r)) / h_i_sq + c_i * w_i * m_i * (1 - a(1, r)) + p_i_c * (w_i * a(1, r) / R(1, r) + c_i * w_i * a(1, r) / f_c - c_i * w_i * (1 - a(1, r)) / f_i_l);
    disp(p(1, r));
    % next iteration, update r
    r = r + 1;
end
lambda_opt = lambda(1, max_iterations);
mu_opt = mu(1, max_iterations);
fprintf("l=%f\n", lambda_opt);
fprintf("m=%f\n", mu_opt);

% compute optimal a_i, t_i
R_i = B/log(2) * lambertw(0, (h_i_sq * (lambda_opt + (1+ mu_opt) * p_i_c)) / (N_o * B * exp(1) * (1+ mu_opt)) - 1/exp(1)) + 1;
a_i = ((T_i * f_i_l) - (c_i * w_i)) / ((w_i * f_i_l / R_i) + (c_i * w_i * f_i_l / f_c) - (c_i * w_i));
t_i = w_i * a_i / R_i;

fprintf("a_i = %f\n", a_i);
fprintf("R_i = %f\n", R_i);
fprintf("t_i = %f\n", t_i);

% compute minimum power consumption
P = t_i * y(N_o, B, d_i / t_i) / h_i_sq + c_i * w_i * m_i * (1 - a_i) + p_i_c * (w_i * a_i / R_i + c_i * w_i * a_i / f_c - c_i * w_i * (1 - a_i) / f_i_l);
fprintf("Power = %f\n", P);

plot(x, p, 'r');
