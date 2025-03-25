clc; clear; close all;

figure; hold on;

function u = steepest_descent(A, f, tol, max_iter)
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
 
    for k=1:max_iter
        r_old = f - A*u_old;
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

u = steepest_descent(A, f, tol, max_iter);

grid on;
tiledlayout(2,1);
tile1=nexttile;
hold(tile1,'on');
%Analytic plot
analytic = 1/pi^2*sin(pi*x)+2/pi*(1-1/pi^2)*x;
%analytic = 1/pi^2*sin(pi*x)+2*x;
plot(tile1, x,analytic, '-r', 'DisplayName', 'Analytic');
legend;
% Plot solution
plot(tile1, x, u, 'b.-', 'DisplayName', 'u(x)');
xlabel('x'); ylabel('u(x)');
title('Solution of Poisson Equation using Steepest Descent Method');
legend;
grid on;

diff = abs(analytic-u);
tile2=nexttile;
plot(tile2,x,diff,'b-', 'DisplayName', 'Error');
title('Error Distribution');
legend;

