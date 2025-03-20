clc; clear; close all;

h_values = [0.1, 0.01, 0.001];

figure; hold on;

for h = h_values

    x = 0:h:1;
    N = length(x) - 2; 

    A = (1/h^2) * (diag(-2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1));
    b = 1 + x(2:end-1);  

    u = A \ reshape(b, [], 1);

    u_full = [0; u; 0];

    if h==0.1 
        plot_color = '-r';
    elseif h==0.01
        plot_color = '-b';
    else
        plot_color = '-og';
    end
    e=reshape(eig(A), [], 1);
    e_full = [0; e ;0];
    plot(x, e_full, plot_color, 'DisplayName', sprintf('h = %.3f', h))
    xlabel('x');
    ylabel('u(x)');
    title(sprintf('Poisons Equation [0,1]: Error Distribution h = 0.1, 0.01, 0.001'));
    legend;
end

