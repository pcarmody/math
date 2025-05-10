function [u, resvec,k] = conjugate_gradient_method(A, f, N, h, tol, max_iter, M)
    u = zeros(N,1);                % Initial guess
    r = f - A*u;                   % Initial residual
    if M == 0
        p = r;                     % Initial direction
        rs_old = r' * r;           % Initial residual norm squared
    else
        M = chol(A,'lower');
        precond = @(r) M' \ (M \ r);
        z = precond(r);
        p=z;                       % preconditioned residual
        rs_old = r' * z;           % Initial residual preconditioned
    end
    res0 = norm(r);                % Initial residual norm
    resvec = zeros(max_iter,1);   % For storing residuals
    
    for k = 1:max_iter
        Ap = A * p;
        alpha = rs_old / (p' * Ap);
        u = u + alpha * p;
        r = r - alpha * Ap;
        res = norm(r) / res0;
        resvec(k) = res;
        
        if res < tol
            fprintf('CG converged at iteration %d with relative residual %.2e\n', k, res);
            break;
        end
        
	if M == 0
            rs_new = r' * r;
            beta = rs_new / rs_old;
            p = r + beta * p;
            rs_old = rs_new;
        else
            z = precond(r);
            rs_new = r' * z;
            beta = rs_new / rs_old;
            p = z + beta * p;
            rs_old = rs_new;
	end
    end
    
end
