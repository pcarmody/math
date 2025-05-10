function x = gauss_seidel(A, b, N, h, tol, max_iter)
    % Solves Ax = b using Gause Segal's iterative method
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
            fprintf('Converged in %d iterations.\n', k);
            conv = conv_loc;
            return;
        end
        x_old = x;  % Update solution
    end
    
    fprintf('Max iterations reached without convergence.\n');
    conv = conv_loc;
end
