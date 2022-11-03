function w = LGL_weights(tau)

N = length(tau)-1;

[~,y] = lepoly(N,tau(2:end-1));

w = [2/(N*(N+1));2./(N*(N+1)*y.^2);2/(N*(N+1))];

end