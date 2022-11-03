%%

function Out = Zermelo_Pseudospectral()

%%
p.ns = 2; % Number of States
p.nu = 1; % Number of Controls

p.t0 = 0; % Initial time
p.tf = 1; % Final time

p.x0 = 0; % Initial condition for y1
p.z0 = 0; % Initial condition for y2

p.nt = 16; % Number of Node Points

%%
p.tau = LGL_nodes(p.nt-1);  % Scaled time horizon
p.D   = LGL_Dmatrix(p.tau); % For defect constraints
p.w   = LGL_weights(p.tau); % For gaussian quadrature

%%
% Discretized variable indices in X = [x,z,u];
p.xi = 1:p.nt;
p.zi = p.nt+1:2*p.nt;
p.ui = 2*p.nt+1:3*p.nt;

%%

X0 = zeros(p.nt*(p.ns+p.nu),1); % Initial Guess (all zeros)

options = optimoptions(@fmincon,'Display','final','MaxFunEvals',1e5);

%%

[sol,val,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x) objfun(x,p),X0,[],[],[],[],[],[],@(x) constraints(x,p),options);

%%

x = sol(p.xi);
z = sol(p.zi);
u = sol(p.ui);

lamx  = lambda.eqnonlin(p.xi)./p.w;
lamz  = lambda.eqnonlin(p.zi)./p.w;

t = (p.tau*(p.tf-p.t0)+(p.tf+p.t0))/2;

Out.t = t;
Out.x = x;
Out.z = z;
Out.u = u;
Out.lamx = lamx;
Out.lamz = lamz;

end

%%
function f = objfun(x,p)
L  = -x(p.xi);
f  = -x(p.xi(end));%(p.tf-p.t0)/2 * dot(p.w,L);
end

%%
function [c,ceq] = constraints(x,p)

xx = x(p.xi);
zz = x(p.zi);
u  = x(p.ui);

%%
c = [u-pi -pi-u];

%%
Y = [xx,zz];
V = 1;
K = 0.3;
F = (p.tf-p.t0)/2*[V*cos(u)+K.*zz,V*sin(u)];
ceq1 = p.D*Y-F;

ceq2 = xx(1)-p.x0;
ceq3 = zz(1)-p.z0;

ceq = [ceq1(:);ceq2;ceq3;];

end

%%
function tau = LGL_nodes(N)

thetak = (4*[1:N]-1)*pi/(4*N+2);
sigmak = -(1-(N-1)/(8*N^3)-(39-28./sin(thetak).^2)/(384*N^4)).*cos(thetak);
    
ze = (sigmak(1:N-1)+sigmak(2:N))/2;
ep = eps*10;

ze1 = ze+ep+1;

while max(abs(ze1-ze))>=ep
    ze1 = ze;
    [dy,y] = lepoly(N,ze);
    ze = ze-(1-ze.*ze).*dy./(2*ze.*dy-N*(N+1)*y);
end

tau = [-1,ze,1]';

end

%%
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

%%

function w = LGL_weights(tau)

N = length(tau)-1;

[~,y] = lepoly(N,tau(2:end-1));

w = [2/(N*(N+1));2./(N*(N+1)*y.^2);2/(N*(N+1))];

end

%%

function [varargout] = lepoly(n,x)

if nargout == 1
   if n == 0, varargout{1} = ones(size(x)); return; end
   if n == 1, varargout{1} = x; return; end
   polylst = ones(size(x)); 
   poly = x;
   
   for k = 2:n
       polyn = ((2*k-1)*x.*poly-(k-1)*polylst)/k;
       polylst = poly;
       poly = polyn;
   end
   varargout{1} = polyn;
end
   
if nargout == 2
    
    if n == 0, varargout{2} = ones(size(x)); varargout{1} = zeros(size(x)); return; end
    if n == 1, varargout{2} = x; varargout{1} = ones(size(x)); return; end
    
    polylst = ones(size(x));
    pderlst = zeros(size(x));
    poly = x;
    pder = ones(size(x));
    
    for k = 2:n
        polyn = ((2*k-1)*x.*poly-(k-1)*polylst)/k;
        pdern = pderlst+(2*k-1)*poly;
        
        polylst = poly;
        poly = polyn;
        
        pderlst = pder;
        pder = pdern;
    end
    varargout{2} = polyn;
    varargout{1} = pdern;
   
end

end
