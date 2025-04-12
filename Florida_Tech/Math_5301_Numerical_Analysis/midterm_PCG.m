% Parameters
N = 128;
L = pi / 2;
h = L / (N + 1);
x = linspace(h, L - h, N)';  % Interior grid points

% Construct RHS vector f
f = sin(pi * x);
f(end) = f(end) + 1 / h^2;  % Adjust for u(pi/2) = 1 boundary

% Construct matrix A (sparse tridiagonal)
e = ones(N,1);
A = spdiags([-e 2*e -e], -1:1, N, N) / h^2;

% Conjugate Gradient Method (manual implementation)
u = zeros(N,1);                % Initial guess
r = f - A*u;                   % Initial residual
p = r;                         % Initial direction
rs_old = r' * r;               % Initial residual norm squared
tol = 1e-4;                    % Relative tolerance
res0 = norm(r);                % Initial residual norm
max_iter = N;
resvec = zeros(max_iter,1);   % For storing residuals

for k = 1:max_iter
    Ap = A * p;
    alpha = rs_old / (p' * Ap);
    u = u + alpha * p;
    r = r - alpha * Ap;
    res = norm(r) / res0;
    resvec(k) = res;
    
    if res < tol
        fprintf('CG converged at iteration %d with relative residual %.2e\n', k, res);
        break;
    end
    
    rs_new = r' * r;
    beta = rs_new / rs_old;
    p = r + beta * p;
    rs_old = rs_new;
end

% Add boundary values
u_full = [0; u; 1];
x_full = linspace(0, L, N + 2)';

% Plot solution
figure;
plot(x_full, u_full, 'b-o')
xlabel('x')
ylabel('u(x)')
title('Solution of -u''''(x) = sin(\pi x) using custom CG method')
grid on

% Plot convergence history
figure;
semilogy(1:k, resvec(1:k), 'r-o')
xlabel('Iteration')
ylabel('Relative Residual')
title('Convergence of Conjugate Gradient Method')
grid on

