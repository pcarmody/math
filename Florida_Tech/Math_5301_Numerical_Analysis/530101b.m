x=linspace(1,4, 37);
y=1./x;
xx=1:0.08:4;
s=spline(x,y);
plot(x,y,s,xx);
