% leapfrog.m
% Leapfrog for u_t + u_x = 0 on [0,1] periodic
% Initial: u0(x) = sin^{80}(pi x)
% dx = 0.01, dt = 0.008 (lambda = 0.8)

clear; close all; clc;

%% Parameters (coarse)
dx = 0.01;
dt = 0.008;
lambda = dt/dx;
L = 1.0;
Nx = round(L/dx);
x = (0:Nx-1)'*dx;        % column vector grid points (0..1-dx)
Tlist = [0.8, 10.0, 20.0];

fprintf('Coarse grid: Nx=%d, dx=%.4f, dt=%.4f, lambda=%.4f\n', Nx, dx, dt, lambda);

%% 1) Von Neumann analysis: dissipation & dispersion
theta = linspace(0, pi, 2000);           % theta = k*dx
a = lambda * sin(theta);

% Solve quadratic G^2 + 2i a G - 1 = 0 -> physical root:
G = sqrt(1 - a.^2) - 1i * a;             % select root with positive cos part
amp = abs(G);                            % |G(theta)|
% numerical frequency: phi = asin(a)  (omega_num * dt = phi)
phi = asin(a);                           % principal branch in [-pi/2, pi/2]
omega_num = phi / dt;                    % numerical angular frequency (rad/time)
omega_exact = (theta / dx);              % exact for ut + ux = 0 is omega = k = theta/dx

% Avoid division by zero at theta=0
ratio = ones(size(omega_num));
idx = abs(omega_exact) > 0;
ratio(idx) = omega_num(idx) ./ omega_exact(idx);
ratio(~idx) = 1;

% Plot dissipation (should be 1)
figure(1); clf;
plot(theta, amp, 'LineWidth', 1.5);
xlabel('\theta = k \Delta x'); ylabel('|G(\theta)|');
title('Leapfrog dissipation curve (|G(\theta)|)'); grid on; xlim([0 pi]);

% Plot dispersion (omega_num / omega_exact)
figure(2); clf;
plot(theta, ratio, 'LineWidth', 1.5);
xlabel('\theta = k \Delta x'); ylabel('\omega_{num}/\omega_{exact}');
title('Leapfrog dispersion: \omega_{num} / \omega_{exact}'); grid on; xlim([0 pi]);

%% Utility functions: Leapfrog with Lax-Wendroff start, Lax-Wendroff full
function u_out = LF_evolve(u0, dx, dt, T)
    % Leapfrog with one Lax-Wendroff start step
    Nx = length(u0);
    lam = dt/dx;
    nsteps = round(T/dt);
    if nsteps == 0
        u_out = u0; return;
    end
    % compute u^1 via Lax-Wendroff
    u = u0;
    up = circshift(u, -1); um = circshift(u, 1);
    u1 = u - 0.5*lam*(up - um) + 0.5*(lam^2)*(up - 2*u + um);
    if nsteps == 1
        u_out = u1; return;
    end
    u_prev = u;       % u^{n-1} = u^0
    u_curr = u1;      % u^{n}   = u^1
    for n = 1:(nsteps-1)
        up = circshift(u_curr, -1);
        um = circshift(u_curr, 1);
        u_next = u_prev - lam*(up - um);   % u^{n+1} = u^{n-1} - lam(...) 
        u_prev = u_curr;
        u_curr = u_next;
    end
    u_out = u_curr;
end

%% Initial condition
u0 = sin(pi*x).^80;

%% Run coarse Leapfrog at requested times
results_coarse = struct();
for k = 1:length(Tlist)
    T = Tlist(k);
    u_num = LF_evolve(u0, dx, dt, T);
    x_shift = mod(x - T, 1);                % analytic shift (u(x,t)=u0(x - t) for ut+ux=0)
    u_exact = sin(pi*x_shift).^80;
    E_err = sqrt(sum((u_num - u_exact).^2) * dx);
    E_exact = sqrt(sum((u_exact.^2)) * dx);
    rel_err = E_err / E_exact;
    results_coarse(k).T = T;
    results_coarse(k).x = x;
    results_coarse(k).u_num = u_num;
    results_coarse(k).u_exact = u_exact;
    results_coarse(k).E_err = E_err;
    results_coarse(k).rel_err = rel_err;
    fprintf('Coarse LF: T=%.2f, L2 err = %.6e, rel = %.6e\n', T, E_err, rel_err);
end

%% Plot comparisons (coarse)
for k = 1:length(Tlist)
    T = Tlist(k);
    figure(10+k); clf;
    plot(results_coarse(k).x, results_coarse(k).u_exact, '-', 'LineWidth', 1.5); hold on;
    plot(results_coarse(k).x, results_coarse(k).u_num, '--', 'LineWidth', 1.2);
    title(sprintf('Leapfrog vs Analytical at t = %.2f (coarse dx=%.3f)', T, dx));
    xlabel('x'); ylabel('u'); legend('Analytical','Leapfrog','Location','Best'); grid on;
    text(0.02, 0.95, sprintf('L2 err = %.3e, rel = %.3e', results_coarse(k).E_err, results_coarse(k).rel_err), ...
        'Units','normalized','VerticalAlignment','top'); hold off;
end

%% PSD of initial condition (coarse and refined) - positive theta side
U0 = fft(u0);
dx_f = 0.005;
dt_f = 0.004;         % keep lambda = 0.8
Nx_f = round(L/dx_f);
x_f = (0:Nx_f-1)'*dx_f;
u0_f = sin(pi*x_f).^80;
psd0 = (abs(U0).^2)/Nx;
kfreqs = (0:Nx-1)'*(2*pi/Nx);      % discrete angular index (mapped to theta ~ k*dx)
npos = floor(Nx/2)+1;
theta_pos = kfreqs(1:npos) * dx;   % approximate theta = k*dx for plotting
figure(60); clf;
semilogy(theta_pos, psd0(1:npos), 'LineWidth', 1.4); hold on;
U0f = fft(u0_f); psdf = (abs(U0f).^2)/Nx_f;
kfreqs_f = (0:Nx_f-1)'*(2*pi/Nx_f);
theta_pos_f = kfreqs_f(1:floor(Nx_f/2)+1) * dx_f;
semilogy(theta_pos_f, psdf(1:floor(Nx_f/2)+1), '--', 'LineWidth', 1.2);
xlabel('\theta = k \Delta x (positive)'); ylabel('PSD (log scale)');
title('Power Spectral Density of initial condition (coarse & refined)');
legend('coarse','refined'); grid on; hold off;
