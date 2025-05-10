% Parameters
N = 128;                % Number of intervals
L = pi / 2;             % Length of domain
h = L / N;              % Grid spacing
x = linspace(0, L, N+1);
x_internal = x(2:N);    % Interior points (excluding boundaries)
n = length(x_internal); % 127

% Right-hand side
f = (h^2) * sin(pi * x_internal)';

% Apply boundary conditions
f(end) = f(end) + 1;    % u(pi/2) = 1, added to the last entry of f

% Create matrix A (tridiagonal: -1, 2, -1)
e = ones(n,1);
A = spdiags([-e 2*e -e], -1:1, n, n);

% Cholesky preconditioner (exact)
Lchol = chol(A, 'lower');    % A = L * L'

% Preconditioning: Solve M^{-1}r using Cholesky
precond = @(r) Lchol' \ (Lchol \ r);

% Initial guess
u = zeros(n, 1);
r = f - A * u;
z = precond(r);              % Preconditioned residual
p = z;
rz_old = r' * z;

% CG Iteration
maxit = 500;
tol = 1e-10;
for k = 1:maxit
    Ap = A * p;
    alpha = rz_old / (p' * Ap);
    u = u + alpha * p;
    r = r - alpha * Ap;
    
    if norm(r) < tol
        fprintf('Converged at iteration %d\n', k);
%        break;
    end
    
    z = precond(r);
    rz_new = r' * z;
    beta = rz_new / rz_old;
    p = z + beta * p;
    rz_old = rz_new;
end

% Add boundary values
u_full = [0; u; 1];

% Plot
plot(x, u_full, 'LineWidth', 2)
xlabel('x'), ylabel('u(x)')
title('1D Poisson Equation Solution using CG + Cholesky Preconditioning')
grid on

