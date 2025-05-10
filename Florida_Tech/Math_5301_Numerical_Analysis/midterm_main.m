clc; clear; close all;

figure; hold on;

%   Define the functions that will converge towards a precise solution

function [x, conv] = gauss_seidel(A, b, N, h, tol, max_iter)
    % Solves Ax = b using Gauss Segal's iterative method
    % Inputs:
    %   A        - Coefficient matrix (NxN)
    %   b        - Right-hand side vector (Nx1)
    %   tol      - Convergence tolerance
    %   max_iter - Maximum number of iterations
    % Output:
    %   x - Solution vector (Nx1)
    
%    N = length(b);              % Number of equations
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
            fprintf('Gauss-Seidel: Converged in %d iterations.\n', k);
            conv = conv_loc;
            return;
        end
        x_old = x;  % Update solution
    end
    
    fprintf('Gauss-Seidel: Max iterations reached without convergence.\n');
    conv = conv_loc;
end

%
% Define the function for CG and CG w/Cholesky Precondition
%    the last parameter, lchol, is a flag.  If set it uses the Cholesky precondition
%

function [u, resvec,k] = conjugate_gradient_method(A, f, N, h, tol, max_iter, lchol)
    u = zeros(N,1);                % Initial guess
    r = f - A*u;                   % Initial residual
    if lchol == 0
        p = r;                     % Initial direction
        rs_old = r' * r;           % Initial residual norm squared
    else
        lchol = chol(A,'lower');%lchol;
        precond = @(r) lchol' \ (lchol \ r);
        z = precond(r);
        p=z;                       % preconditioned residual
        rs_old = r' * z;           % Initial residual preconditioned
    end
    res0 = norm(r);                % Initial residual norm
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
        
	if lchol == 0
            rs_new = r' * r;
            beta = rs_new / rs_old;
            p = r + beta * p;
            rs_old = rs_new;
        else
            z = precond(r);%lchol' \ (lchol \ r)%precond(r);
            rs_new = r' * z;
            beta = rs_new / rs_old;
            p = z + beta * p;
            rs_old = rs_new;
	end
    end
    
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

% Solve using Gauss-Segal method
tol = 1e-4;  % Convergence tolerance
max_iter = 1000;  % Maximum iterations

[gs_plot, gs_conv] = gauss_seidel(A, f, N, h, tol, max_iter);

[cg_plot, cg_resvec, k] = conjugate_gradient_method(A, f, N, h, tol, max_iter, 0);

[cgc_plot, cgc_resvec, kc] = conjugate_gradient_method(A, f, N, h, tol, max_iter, 1);
x_full = linspace(0, pi/2, N )';

%Analytic plot
analytic = 1/pi^2*sin(pi*x)+2/pi*(1-1/pi^2)*x;
grid on;
tiledlayout(2,1);
tile1=nexttile;
hold(tile1,'on');

horizontal = 1:length(gs_plot);
plot(tile1, x_full, analytic, '-g', 'DisplayName', 'Analytic');
plot(tile1, x_full, gs_plot, '-r', 'DisplayName', 'Gauss-Segal');
plot(tile1, x_full, cg_plot, 'b o', 'DisplayName', 'Conjugate Gradient');
plot(tile1, x_full, cgc_plot, 'k', 'DisplayName', 'Conjugate Gradient w/Preconditioning');
xlabel('x'); ylabel('u(x)');
title('Poisson Equation using Gauss-Segal Method');
legend;
grid on;

tile2=nexttile;
hold(tile2,'on');
gs_diff = abs(analytic-gs_plot);
cg_diff = abs(analytic-cg_plot);
cgc_diff = abs(analytic-cgc_plot);
plot(tile2, x_full, gs_diff, '-r', 'DisplayName', 'Analytic - GS');
plot(tile2, x_full, cg_diff, '-b', 'DisplayName', 'Analytic - CG');
plot(tile2, x_full, cgc_diff, '-k', 'DisplayName', 'Analytic - CGw/P');
title('Error Distribution');
legend;

figure;
semilogy(1:k-1, cg_resvec(1:k-1), 'b-o', 'DisplayName', 'Convergence of CG');
xlabel('Iteration')
ylabel('Relative Residual')
title('Convergence of Conjugate Gradient Method')
legend;
grid on
