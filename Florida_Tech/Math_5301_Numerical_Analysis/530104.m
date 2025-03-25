clc; clear; close all;

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
plot(x, u_full, '-b', 'DisplayName', 'Central Difference');
xlabel('x');
ylabel('u(x)');
title(sprintf('Poisons Equation: Solution using mesh=128'));
legend;
grid on;

