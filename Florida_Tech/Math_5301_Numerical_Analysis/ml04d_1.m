clc; clear; close all;

figure; hold on;

function [conv, cnt] = steepest_descent(A, f, tol, max_iter)
    % Solves Ax = b using steepest descent method
    % Inputs:
    %   A        - Coefficient matrix (NxN)
    %   b        - Right-hand side vector (Nx1)
    %   tol      - Convergence tolerance
    %   max_iter - Maximum number of iterations
    % Output:
    %   x - Solution vector (Nx1)
    
    N = length(f);              % Number of equations
    u = zeros(N, 1);            % Initial guess (zero vector)
    u_old = u;
    conv = 1:10000;
 
    for k=1:max_iter
        r_old = f - A*u_old;
        conv(k) = norm(r_old);
        cnt = k;
        
        if norm(r_old) < tol
            fprintf('Converged in %d iterations.\n', k);
            return;
        end
        alpha_old = (r_old.' * r_old)/(r_old.' * A * r_old);
        u = u_old +  alpha_old*r_old;
        u_old = u;
    end
    
    fprintf('Max iterations reached without convergence.\n');
end
% Define problem parameters
N = 128;                  % Number of internal grid points
h = (pi/2) / (N+1);       % Grid spacing
x = linspace(h, pi/2-h, N)';  % Grid points
f = h^2* sin(pi * x);    % Right-hand side vector

% Construct tridiagonal matrix A
A = 2 * eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);

% Modify last element of f to include boundary condition u(pi/2) = 1
f(1) = 0; f(N) = 1;

% Solve using steepest descent method
tol = 1e-4;  % Convergence tolerance
max_iter = 100000;  % Maximum iterations

[conv, size] = steepest_descent(A, f, tol, max_iter);

output = conv(1:size);
grid on;
yscale log;
horizontal = 1:length(output);%linspace(1, cnt, 1);
plot(horizontal, output, '-k', 'DisplayName', 'Rate of Convergence');
xlabel('x'); ylabel('tolerance');
title('Convergenace of Poisson Equation using Stepest Descent Method');
legend;

