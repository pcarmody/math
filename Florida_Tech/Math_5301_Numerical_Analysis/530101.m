x=linspace(0,pi/2,100);
sinx=sin(x);
p2=-0.3357*x.^2+0.9003*x;
error=(1/6)*x.*(x-pi/4).*(x-pi/2);
errhi=p2+error;
errlo=p2-error;
plot(x,sinx);
plot(x,p2);
plot(x,errhi);
plot(x,errlo);
legend("sin(x)", "P2", "Error Hi", "Error Low");
