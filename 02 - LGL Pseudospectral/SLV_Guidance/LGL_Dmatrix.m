function D = LGL_Dmatrix(tau)

N = length(tau)-1;

n = N+1;

if n == 0, D = []; return; end

xx = tau;
y  = lepoly(n-1,xx);

D = (xx./y)*y'-(1./y)*(xx.*y)';

D = D+eye(n);
D = 1./D;
D = D-eye(n);

D(1,1) = -n*(n-1)/4;
D(n,n) = -D(1,1);

end