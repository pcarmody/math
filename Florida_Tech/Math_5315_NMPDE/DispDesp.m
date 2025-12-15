% DispDesp.m
% FTFS for u_t - u_x = 0 on [0,1] with periodic BC
% Delta x = 0.01, Delta t = 0.008 (lambda = dt/dx = 0.8)
% Produces: dissipation & dispersion curves, numerical vs analytical solutions,

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
