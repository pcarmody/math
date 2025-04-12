clc; clear; close all;

figure; hold on;

function x = gause_segal(A, b, h, tol, max_iter)
    % Solves Ax = b using Gause Segal's iterative method
    % Inputs:
    %   A        - Coefficient matrix (NxN)
    %   b        - Right-hand side vector (Nx1)
    %   tol      - Convergence tolerance
    %   max_iter - Maximum number of iterations
    % Output:
    %   x - Solution vector (Nx1)
    
    N = length(b);              % Number of equations
    x = b;%zeros(N, 1);            % Initial guess (zero vector)
    x_old = x;                  % Store previous iteration
    h2=h^2;
    conv_loc = 1:10000;
    conv=conv_loc;
    cnt=0;
 
    for k=1:max_iter
        for i=2:N-1
%            x(i) = 1/2*( x_old(i-1) + x_old(i+1) - h2*b(i) ); % Jacobi
            x(i) = 1/2*( x(i-1) + x_old(i+1) - h2*b(i) );  % Gaus-Segal
        end
        
        % Check for convergence
%        if norm(x - x_old, 2) < tol
         conv_loc(k) = norm(x-x_old, Inf);
         cnt = k;
        
        if norm(x - x_old, Inf) < tol
            fprintf('Converged in %d iterations.\n', k);
            conv = conv_loc;
            return;
        end
        x_old = x;  % Update solution
    end
    
    fprintf('Max iterations reached without convergence.\n');
    conv = conv_loc;
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

% Solve using Gause-Segal method
tol = 1e-4;  % Convergence tolerance
max_iter = 1000;  % Maximum iterations

u = gause_segal(A, f, h, tol, max_iter);

output = u;%u(1:cnt);
grid on;
%yscale log;
horizontal = 1:length(output);%linspace(1, cnt, 1);
plot(horizontal, output, '-k', 'DisplayName', 'Gause-Segal');
xlabel('x'); ylabel('u(x)');
title('Convergenace of Poisson Equation using Gause-Segal Method');
legend;
grid on;

