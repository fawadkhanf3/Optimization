%% Zermelo's Problem (With Collocation Method)

function Out = Zermelo_Collocation()

%%
global V tf K xdot zdot N h

N  = 40; % Number of Segments
K  = 0.3;
V  = 1;
tf = 1;

xdot = @(x,z,u) V*cos(u) + K*z;
zdot = @(x,z,u) V*sin(u);

%% Simulation Setup

a = 0;
b = tf;
h = (b-a)/N;

t = linspace(a,b,N+1)';

% Initial Guess
x  = ones(N+1,1);
z  = ones(N+1,1);
u  = zeros(N+1,1);
ud = zeros(N+1,1);
C2 = ones(N,1);
C3 = ones(N,1);

xx0 = [x ;z; u ; ud ; C2; C3];

cost_func = @(xx) -xx(N+1); % Maximizing final value of x

opt     = optimset('fmincon');
options = optimset(opt,'TolX',0.0001,'TolCon',0.0001,...
    'MaxIter',5000,'MaxFunEvals',1e6,...
    'Display','Final');

%% Solution

[xx,val,exitflag,output,lambda,grad,hessian] = ...
    fmincon(cost_func,xx0,[],[],[],[],[],[],@constraints,options);

%% Results

n     = N+1;
x     = xx(1:n);
z     = xx(n+1:2*n);
gamma = xx(2*n+1:3*n);
lamx  = lambda.eqnonlin(1:N);
lamz  = lambda.eqnonlin(N+1:2*N);

Out.t = t;
Out.x = x;
Out.z = z;
Out.u = gamma;
Out.lamx = lamx;
Out.lamz = lamz;

end

%% Constraints

function [c,d] = constraints(xx)

global N h xdot zdot

% Assigning States; at all node points (N+1)
n  = N+1;
x  = xx(1:n);
z  = xx(n+1:2*n);
u  = xx(2*n+1:3*n);
ud = xx(3*n+1:4*n);

%% Inequality Constraints on u
c = [u-pi -pi-u]; % -pi <= u <= pi

%% Equality Constraints

% For N-Segments
C0 = u(1:N);
C1 = ud(1:N);
C2 = xx(4*n+1:4*n+N);
C3 = xx(4*n+N+1:4*n+N+N);

% RHS of ODE at all node points (N+1)
f1 = xdot(x,z,u);
f2 = zdot(x,z,u);

% States & Derivatives at Midpoints (N)
xMP  = zeros(N,1); zMP  = zeros(N,1);
xdMP = zeros(N,1); zdMP = zeros(N,1);
for i = 1:N
    % States at Midpoint
    xMP(i) = (x(i)+x(i+1))/2 + h*(f1(i) - f1(i+1))/8;
    zMP(i) = (z(i)+z(i+1))/2 + h*(f2(i) - f2(i+1))/8;
    
    % Derivatives at Midpoint
    xdMP(i) = -3*(x(i)-x(i+1))/2/h -(f1(i) + f1(i+1))/4;
    zdMP(i) = -3*(z(i)-z(i+1))/2/h -(f2(i) + f2(i+1))/4;
end

% Control at MidPoints (S = 1/2) (N)
S   = 1/2;
uMP = C0 + C1.*S + C2.*S^2 + C3.*S^3;

% R.H.S of ODE at Midpoints (N)
f1MP = xdot(xMP,zMP,uMP);
f2MP = zdot(xMP,zMP,uMP);

% Defect Equations (Derivative at M.P == R.H.S of ODE evaluated at M.P)
d(1:N)     = f1MP - xdMP;
d(N+1:2*N) = f2MP - zdMP;

% Continuity Constraints (S=1 of Previous Node == S=0 of Next Node)
for i=1:N
    d(2*N+i) = (C0(i)+C1(i)+C2(i)+C3(i)) - u(i+1);
    d(3*N+i) = (C1(i)+2*C2(i)+3*C3(i))   - ud(i+1);
end

for i=1:N-1
    d(4*N+i) = (2*C2(i)+6*C3(i)) - 2*C2(i+1);
end

% Initial Value Constraint
d(5*N)   = x(1) - 0;
d(5*N+1) = z(1) - 0;

end