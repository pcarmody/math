clc; clear; close all;

figure; hold on;

function x = jacobi_method(A, b, tol, max_iter)
    % Solves Ax = b using Jacobi's iterative method
    % Inputs:
    %   A        - Coefficient matrix (NxN)
    %   b        - Right-hand side vector (Nx1)
    %   tol      - Convergence tolerance
    %   max_iter - Maximum number of iterations
    % Output:
    %   x - Solution vector (Nx1)
    
    N = length(b);              % Number of equations
    x = zeros(N, 1);            % Initial guess (zero vector)
    x_old = x;                  % Store previous iteration
    D = diag(A);                % Diagonal elements of A
    R = A - diag(D);            % Remaining part of A (L+U)
    
    for k = 1:max_iter
        x = (b - R * x_old) ./ D;  % Jacobi iteration formula
        
        % Check for convergence
        if norm(x - x_old, inf) < tol
            fprintf('Converged in %d iterations.\n', k);
            return;
        end
        x_old = x;  % Update solution
    end
    
    fprintf('Max iterations reached without convergence.\n');
end
% Define problem parameters
N = 128;                  % Number of internal grid points
h = (pi/2) / (N+1);       % Grid spacing
x = linspace(h, pi/2-h, N)';  % Grid points
f = h^2 * sin(pi * x);    % Right-hand side vector

% Construct tridiagonal matrix A
A = 2 * eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);

% Modify last element of f to include boundary condition u(π/2) = 1
f(N) = f(N) + 1;

% Solve using Jacobi method
tol = 1e-4;  % Convergence tolerance
max_iter = 10000;  % Maximum iterations

u = jacobi_method(A, f, tol, max_iter);

% Plot solution
plot(x, u, 'b.-');
xlabel('x'); ylabel('u(x)');
title('Solution of Poisson Equation using Jacobi Method');
grid on;

