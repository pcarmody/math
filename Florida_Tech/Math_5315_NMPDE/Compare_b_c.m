% Compare_b_c.m

clear; close all; clc;

%% Parameters
dx = 0.01;
dt = 0.008;
lambda = dt/dx;            % Courant number
L = 1;
Nx = round(L/dx);
x = (0:Nx-1)'*dx;          % column vector
T = 1.0;
nsteps = round(T/dt);      % should be integer (125)
if abs(nsteps*dt - T) > 1e-12
    error('T not an integer multiple of dt');
end

%% Helper: FTFS step function (periodic)
ftfs_step = @(u) (1-lambda)*u + lambda * circshift(u,-1);  % forward shift is j+1 -> circshift -1

%% ------ CASE (b): u0 = sin^{40}(pi x) ------
u0_b = sin(pi*x).^40;
u = u0_b;
for n = 1:nsteps
    u = ftfs_step(u);
end
u_num_b = u;
% Analytical solution: shift by +t (because u(x,t) = u0(x + t) for u_t - u_x = 0)
x_shift = mod(x + T, 1);
u_exact_b = sin(pi*x_shift).^40;

% Energy (L2) (using dx Riemann sum)
E0_b = sum((u0_b).^2) * dx;
En_b = sum((u_num_b).^2) * dx;
frac_diss_b = 1 - En_b / E0_b;

%% Plot analytical vs numerical for (b)
figure(3); clf;
plot(x, u_exact_b, '-', 'LineWidth', 1.6, 'DisplayName', 'Analytical (t=1)');
hold on;
plot(x, u_num_b, '--', 'LineWidth', 1.6, 'DisplayName', 'Numerical FTFS (t=1)');
xlabel('x'); ylabel('u(x,1)');
title('Case (b): u_0(x) = sin^{40}(\pi x)  -- Analytical vs FTFS (t=1)');
legend('Location','Best');
grid on;

fprintf('\nCase (b) (sin^{40}(pi x)):\n  E0 = %.6e, E_num(t=1) = %.6e, fractional dissipation = %.6e\n',...
    E0_b, En_b, frac_diss_b);

%% ------ CASE (c): square pulse on [0.4, 0.6] ------
u0_c = zeros(Nx,1);
% careful: include index where x in [0.4,0.6] inclusive
idx = (x >= 0.4) & (x <= 0.6);
u0_c(idx) = 1.0;

u = u0_c;
for n = 1:nsteps
    u = ftfs_step(u);
end
u_num_c = u;
% Analytical shift
u_exact_c = double((x_shift >= 0.4) & (x_shift <= 0.6));

% Energy
E0_c = sum((u0_c).^2) * dx;
En_c = sum((u_num_c).^2) * dx;
frac_diss_c = 1 - En_c / E0_c;

%% Plot analytical vs numerical for (c)
figure(4); clf;
plot(x, u_exact_c, '-', 'LineWidth', 1.6, 'DisplayName', 'Analytical (t=1)');
hold on;
plot(x, u_num_c, '--', 'LineWidth', 1.6, 'DisplayName', 'Numerical FTFS (t=1)');
xlabel('x'); ylabel('u(x,1)');
title('Case (c): Square pulse on [0.4,0.6]  -- Analytical vs FTFS (t=1)');
legend('Location','Best');
grid on;

fprintf('\nCase (c) (square pulse):\n  E0 = %.6e, E_num(t=1) = %.6e, fractional dissipation = %.6e\n',...
    E0_c, En_c, frac_diss_c);

