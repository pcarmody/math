clc; clear; close all;

h_values = [0.1, 0.01, 0.001];

figure; hold on;

function dsc = centered_difference_scheme(mesh)

    h = 1/mesh;
    x = 0:h:1;
    N = length(x) - 2; 

    % Construct finite difference matrix A
    A = (1/h^2) * (diag(-2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1));
    b = 1 + x(2:end-1);  

    % Solve the linear system A*u = b
    u = A \ reshape(b, [], 1);

    % Include boundary values (assuming Dirichlet B. C. u(0) = u(1) = 0)
    dsc = [0; u; 0];

end

u_full= centered_difference_scheme(128);
% Plot the solution
plot(x, u_full, '-b', 'DisplayName', sprintf('h = %.3f', h))
xlabel('x');
ylabel('u(x)');
title(sprintf('Poisons Equation [0,1]: Solution using h = 0.1, 0.01, 0.001'));
legend;

