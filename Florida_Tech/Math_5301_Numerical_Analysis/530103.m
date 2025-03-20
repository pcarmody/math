clc; clear; close all;

h_values = [0.1, 0.01, 0.001];

figure; hold on;

for h = h_values

    x = 0:h:1;
    N = length(x) - 2; 

    % Construct finite difference matrix A
    A = (1/h^2) * (diag(-2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1));
    b = 1 + x(2:end-1);  

    % Solve the linear system A*u = b
    u = A \ reshape(b, [], 1);

    % Include boundary values (assuming Dirichlet B. C. u(0) = u(1) = 0)
    u_full = [0; u; 0];

    % Plot the solution
    if h==0.1 
        plot_color = '-or';
    elseif h==0.01
        plot_color = '-b';
    else
        plot_color = '-g';
    end
    plot(x, u_full, plot_color, 'DisplayName', sprintf('h = %.3f', h))
    xlabel('x');
    ylabel('u(x)');
    title(sprintf('Poisons Equation [0,1]: Solution using h = 0.1, 0.01, 0.001'));
    legend;
end

