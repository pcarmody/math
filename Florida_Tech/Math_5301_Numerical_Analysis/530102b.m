n=6;
c = 1:1:n;
xc = -1.5*cos(((2.*c-1)/10)*pi); % our Chebeschev Nodes
%xc(3)=0;
fx = abs(xc);
hold on;
plot(xc,fx);

p=polyfit(xc,fx,4);
x=linspace(-1.5,1.5,100);
plot(x,polyval(p,x),'--');
legend("Function", "Polynomial");
