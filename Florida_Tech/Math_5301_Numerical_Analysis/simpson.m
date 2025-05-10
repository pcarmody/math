
global pdata inc;
pdata = createArray(1,100);

function [new_est, depth, iterat] = simpson_third(est, a, b, level, itr, tol)

	% est: estimate of current iteration
	% a: left boundary
	% b: right boundary
	% level: how depend has recursion gone
	% itr: how many iterations have been executed
	% tol:  tolerance level

	% set default responses
	global pdata inc;

	depth = level;
	iterat = itr+1;
	new_est = est;
	if a==b		% we have no interval to consider
		return;
	end

	m = floor((b-a)/2);	% midpoint

	if m <= 1	% cannot recurse any deeper
		return;
	end
	m=m+a;

	% find the new estimate

	new_est = ((b-a)*inc)/6*( pdata(a) + 4 * pdata(m) + pdata(b) );

	if abs(est - new_est) <= tol % is this estimate close enough?
		return;
	end

	new_est = new_est/2.0;
	[left, depth_left, left_itr] = simpson_third(new_est, a, m, level+1, itr, tol);
	[right, depth_right, right_itr] = simpson_third(new_est, m, b, level+1, itr, tol);
	new_est = left+right;
	depth = max(depth_left, depth_right);
	iterat = left_itr + right_itr;

end

n=1000;

% Test 1: y=x, [0, 4]
inc = 4/n;
for i = 1:n
	pdata(i) = 1;%(inc*i)^2;
end;

[t1_est, t1_dept, t1_rtr] = simpson_third(0, 1, n, 0, 0, inc);
fprintf('Test 1: estimate = %2.3f, Depth=%d, Iterations=%d\n', t1_est, t1_dept, t1_rtr);

% Test 2: y=x^2, [0, 4]
inc = 4/n;
for i = 1:n
	pdata(i) = (inc*i)^2;
end;

[t2_est, t2_dept, t2_rtr] = simpson_third(22, 1, n, 0, 0, inc);
fprintf('Test 2: estimate = %2.3f, Depth=%d, Iterations=%d\n', t2_est, t2_dept, t2_rtr);

% Test 3: y=|sin(x)|, [0, pi/2]
inc = pi/2/n;
for i = 1:n
	pdata(i) = abs(sin(i*inc));
end;

[t3_est, t3_dept, t3_rtr] = simpson_third(22, 1, n, 0, 0, inc);
fprintf('Test 3: estimate = %2.3f, Depth=%d, Iterations=%d\n', t3_est, t3_dept, t3_rtr);

% Test 4: y=sin(5x)+2, [0, pi/2]
inc = pi/2/n;
for i = 1:n
	pdata(i) = sin(11*i*inc)+1;
end;

[t4_est, t4_dept, t4_rtr] = simpson_third(22, 1, n, 0, 0, inc);
fprintf('Test 4: estimate = %2.3f, Depth=%d, Iterations=%d\n', t4_est, t4_dept, t4_rtr);

% test 5: complex partitions
inc = pi/2/n;
for i = 1:n/2;
	pdata(i) = inc*i;
end
for i=n/2:n;
	pdata(i) = sin(11*i*inc) + 1;
end

[t5_est, t5_dept, t5_rtr] = simpson_third(22, 1, n, 0, 0, inc);
fprintf('Test 5: estimate = %2.3f, Depth=%d, Iterations=%d\n', t5_est, t5_dept, t5_rtr);
%[est, dept, itr] = simpson_third(pi/2, 1, n, 1, 1, 1e-3)
