p.t0 = 0; % Initial time
p.tf = 50; % Final time

p.r0     = 1.1;
p.theta0 = 0;
p.u0     = 0;
p.v0     = 1/sqrt(1.1);

p.A  = 0.01;

p.ns = 4; % Number of States
p.nu = 1; % Number of Controls

p.nt = 64; % Number of Node Points

% Discretized variable indices in X = [x,z,u];
p.ri     = 1:p.nt;
p.thetai = p.nt+1:2*p.nt;
p.ui     = 2*p.nt+1:3*p.nt;
p.vi     = 3*p.nt+1:4*p.nt;
p.Ei     = 4*p.nt+1:5*p.nt;

%% LGL Nodes, Differential Matrix and Weights

p.tau = LGL_nodes(p.nt-1);  % Scaled time horizon; tau = [-1,+1]
p.D   = LGL_Dmatrix(p.tau); % For defect constraints
p.w   = LGL_weights(p.tau); % For gaussian quadrature

%%


load Solution.mat

t = (p.tau*(p.tf-p.t0)+(p.tf+p.t0))/2;

r     = sol(p.ri);
theta = sol(p.thetai);
u     = sol(p.ui);
v     = sol(p.vi);
E     = sol(p.Ei);

plot(t,E,'r.-');grid;shg

mm = [t E];

x0 = [p.r0;p.theta0;p.u0;p.v0];
tspan = [p.t0 p.tf];


options = odeset('MaxStep',0.5);
[t,y] = ode45(@(t,x) dots(t,x,mm,p),tspan,x0,options);

r = y(:,1);
theta = y(:,2);
u = y(:,3);
v = y(:,4);

figure(1);hold on;grid on
plot(r.*cos(theta),r.*sin(theta),'b.-');

figure(2);hold on;grid on
plot(t,u,'ko-',t,v,'go-');
xlabel('t [sec]');
ylabel('Optimal Velocities');
legend('u','v');