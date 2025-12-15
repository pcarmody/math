% ftfs_wave_analysis.m
% FTFS for u_t - u_x = 0 on [0,1] with periodic BC
% Delta x = 0.01, Delta t = 0.008 (lambda = dt/dx = 0.8)
% Produces: dissipation & dispersion curves, numerical vs analytical solutions,
% PSD of initial conditions, and prints L2-energy dissipation.

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

fprintf('Nx = %d, dx = %.4f, dt = %.4f, lambda = %.4f, nsteps = %d\n', Nx, dx, dt, lambda, nsteps);

%% FTFS operator (for Fourier analysis)
theta = linspace(0, pi, 2000);   % theta = k*dx (0..pi)
G = 1 - lambda + lambda*exp(1i*theta);  % amplification per time-step
amp = abs(G);                       % dissipation per step
phase = angle(G);                   % numerical phase per step
omega_num = -phase / dt;            % numerical angular freq (rad / unit time)
omega_exact = -theta / dx;          % exact angular freq (for wave u_t - u_x = 0 -> omega = -k)
ratio = ones(size(omega_num));
nz = abs(omega_exact) > 1e-16;
ratio(nz) = omega_num(nz) ./ omega_exact(nz);

%% Plot dissipation |G(theta)| vs theta
figure(1); clf;
plot(theta, amp, 'LineWidth', 1.6);
xlabel('\theta = k \Delta x');
ylabel('|G(\theta)| (amplification magnitude per step)');
title('FTFS dissipation curve: |G(\theta)| vs \theta');
grid on;
xlim([0 pi]);

%% Plot dispersion: omega_num/omega_exact vs theta
figure(2); clf;
plot(theta, ratio, 'LineWidth', 1.6);
xlabel('\theta = k \Delta x');
ylabel('\omega_{num} / \omega_{exact}');
title('FTFS dispersion: numerical/exact frequency ratio');
grid on;
xlim([0 pi]);

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

%% ----- PSD analysis of initial conditions -----
% Compute FFTs and PSD (one-sided positive theta 0..pi)
U_b = fft(u0_b);
U_c = fft(u0_c);
psd_b = (abs(U_b).^2) / Nx;
psd_c = (abs(U_c).^2) / Nx;

% Frequency indices for FFT: use fftfreq-like mapping
k = (0:Nx-1)';                          % discrete indices
% physical angular wavenumber k_phys = 2*pi*(k/Nx)/dx? Simpler: theta = k*2*pi/Nx * dx obviously maps to [-pi,pi)
% Use fftshift to center; but we'll just reorder to ascending theta from 0..2pi and map to -pi..pi.
theta_all = (2*pi/Nx) * (0:Nx-1)';     % angular wavenumber * dx? This is in radians per grid cell
% But theta = k*dx in our notation; map k_idx->theta_k = (k if k<=Nx/2 else k-Nx)*2*pi/Nx * dx
k_shift = ifftshift(0:Nx-1) - floor(Nx/2);  % not used directly in plotting; we'll instead produce positive half
% We'll plot positive-theta modes 0..pi (first Nx/2+1 entries)
npos = floor(Nx/2) + 1;
% Build theta positive
theta_pos = (0:(npos-1)) * (2*pi/Nx);   % 0 .. pi approximately
psd_b_pos = psd_b(1:npos);
psd_c_pos = psd_c(1:npos);

figure(5); clf;
semilogy(theta_pos, psd_b_pos, '-', 'LineWidth', 1.6, 'DisplayName', 'PSD sin^{40}(\pi x)');
hold on;
semilogy(theta_pos, psd_c_pos, '--', 'LineWidth', 1.6, 'DisplayName', 'PSD square pulse');
xlabel('\theta = k \Delta x  (0 to \pi)');
ylabel('Power spectral density (log scale)');
title('PSD of initial conditions (positive \theta modes)');
legend('Location','Best');
grid on;

%% Also plot net damping factor after nsteps: |G(theta)|^{nsteps}
G_full = (1 - lambda + lambda*exp(1i*theta)).^nsteps;
amp_full = abs(G_full);

figure(6); clf;
plot(theta, amp_full, 'LineWidth', 1.6);
xlabel('\theta = k \Delta x');
ylabel(['|G(\theta)|^{', num2str(nsteps), '}  (net amplification to t=1)']);
title('Net per-mode amplification (damping) after t=1 (FTFS)');
grid on;
xlim([0 pi]);

%% Print summary
fprintf('\nSummary:\n  lambda = %.3f, Nx = %d, nsteps = %d\n', lambda, Nx, nsteps);
fprintf('  Case (b) fractional dissipation = %.6e\n', frac_diss_b);
fprintf('  Case (c) fractional dissipation = %.6e\n', frac_diss_c);

