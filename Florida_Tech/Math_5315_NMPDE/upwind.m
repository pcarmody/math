
% switch on frame saving
saveframes=0;

% solve the advection equation: du/dt + du/dx = 0.0
% using first order upwinding 

% number of data points
M = 100;

% periodic interval dimensions (i.e. x \in [a,b))
a = 0; b =  1;

% data point spacing
% dx = (b-a)/M;
dx = 0.01;
M = (b-a)/dx;  % reset the number of data points

% coordinates of data points ("nodes")
x = linspace(a, b-dx, M);
u = zeros(1,M);

% start time
StartTime = 0;
EndTime = 2;

% dt chosen small to make space errors dominant
% dt = 0.55*dx;
dt = 0.008;

% compute number of time steps
Nsteps = ceil((EndTime-StartTime)/dt);

% modify dt to end exactly at EndTime
dt = (EndTime-StartTime)/Nsteps;

% initial condition 
for i=1:1:M
  u(i) = exactsolution(x(i));
end
t = StartTime; % set time variable
fram = 1;
for n=1:Nsteps % begin time stepping
    
    um1 = [u(M),u(1:M-1)];       % u_{m-1}
    up1 = [u(2:M),u(1)];
    u1 = [u(2:M),u(1)];
    
     % update solution 
%   unew = u-0.5*dt/dx*(up1-um1)+0.5*dt*dt*(up1+um1-2*u)/(dx*dx); % Wave Equation
%   unew = u-dt/dx*(u-um1);
   unew = u-dt/dx*(up1-u);  %  Transport Equation
   unew = u-dt/dx*(u-um1);
       
    u = unew; % finish one time step
    t  = t+dt;   % update time
    
    if(mod(n,20)==0 |n==Nsteps) % selective plotting
    
        plot(x, u, 'b-o', 'LineWidth', 1.2,'MarkerSize', 6); % plot numerical solution
        hold on;
        xplot = linspace(a, b, 1000); % plot exact solution at lots of points
        xmod = mod(xplot-t,b); % coordinate for exact pulse
        uexact = zeros(1,1000);
        for i=1:1:1000
          uexact(i) = exactsolution(xmod(i));
        end
        plot(xplot, uexact, 'k-', 'LineWidth', 2);hold off;
        %axis([a b -.5 1.5]);
        legend(sprintf('Numerical solution (t=%g)',t), 'Exact solution');
        set(gca,'fontsize',18);

        drawnow; pause(0.05);
        if(saveframes == 1)        % save frame to file
            if(fram<10)            fname = sprintf('anim_000%d.ppm', fram);
            elseif(fram<100)   fname = sprintf('anim_00%d.ppm', fram);
            elseif(fram<1000) fname = sprintf('anim_0%d.ppm', fram);
            end
            print('-dppm', fname);
            fram = fram+1;
        end
    end
end

clear
clc
