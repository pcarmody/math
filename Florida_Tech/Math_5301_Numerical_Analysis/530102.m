x = linspace(-3,3,60); 
y = 1./(1+5*x.^2);

% initialize the working arrays
alpha = linspace(-3,3,60);
l = linspace(-3,3,60);
u = linspace(-3,3,60);
z = linspace(-3,3,60);

% initialize the return values

a = y;
b = linspace(-3,3,60);
c = linspace(-3,3,60);
d = linspace(-3,3,60);

h = 0.1;
n = 61;
for i = 2:1:n-2
	alpha(i) = 30*( y(i+1) - 2*y(i) + y(i-1) );
end
l(1) = 1;
u(1) = 0;
z(1) = 0;

for i = 2:1:n-1
	l(i) = 0.4 - 0.1*u(i-1);
	u(i) = 0.1/l(i);
	z(i) = (alpha(i) - 0.1*z(i-1)) / l(i);
end

l(n)=1;
z(n)=0;
c(n)=0;

% plot f(x) in blue
plot(x,y,'-b');
hold on;

for j = n-2:-1:1
	c(j) = z(j) - u(j)*c(j+1);
	b(j) = ( a(j+1) - a(j) )/0.1 - 0.1*(c(j+1) + 2*c(j) )/3;
	d(j) = ( c(j+1) - c(j) )/30;
end

%now plot each of the S_j splines within their intervals

for j=2:1:n-2
	left = -3 + (j-1)*0.1;
	right = -3 + j*0.1;
	xj=[left:0.01:right];
	p = [a(j) b(j) c(j) d(j)];
	plot(xj, polyval(p, xj));
	%yj=a(j) + b(j).*xj + c(j)*xj.^2 + d(j)*xj.^3;
	%plot(xj,yj,"-r");
	clear xj yj p;
end

