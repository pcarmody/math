clc; clear; close all;

h_values = [0.1, 0.01, 0.001];

figure; hold on;

function [x, dsc] = centered_difference_scheme(mesh)

    h = 1/mesh;
    x = 0:h:pi/2;
    N = length(x) - 2; 

    % Construct finite difference matrix A
    A = (1/h^2) * (diag(-2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1));
    b = sin(pi*x(2:end-1));  

    % Solve the linear system A*u = b
    u = A \ reshape(b, [], 1);

    % Include boundary values (assuming Dirichlet B. C. u(0) = u(1) = 0)
    dsc = [0; u; 0];

end

[x, u_full]= centered_difference_scheme(128);
analytic = (-1/pi^2)*sin(pi*x);
% Plot the solution
plot(x, u_full, '-b', 'DisplayName', sprintf('h = %.3f', 1/128))
plot(x, analytic, '-r', 'DisplayName', 'analytic')
xlabel('x');
ylabel('u(x)');
title(sprintf('Poisons Equation [0,1]: Solution using h = 0.1, 0.01, 0.001'));
legend;

