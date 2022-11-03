%% Zermelo Problem (Solved By LGL Pseudospectral Method)

clear;close all;clc
format long g
set(0,'DefaultLineLineWidth',2);

%% Problem Parameters

p.ns = 2; % Number of States
p.nu = 1; % Number of Controls

p.t0 = 0; % Initial time
p.tf = 1; % Final time

p.x0 = 0; % Initial condition for y1
p.z0 = 0; % Initial condition for y2

p.V  = 1; % Constant Velocity (1 m/s)
p.K  = 0.3;

p.nt = 16; % Number of Node Points

% Discretized variable indices in X = [x,z,u];
p.xi = 1:p.nt;
p.zi = p.nt+1:2*p.nt;
p.ui = 2*p.nt+1:3*p.nt;

%% LGL Nodes, Differential Matrix and Weights

p.tau = LGL_nodes(p.nt-1);  % Scaled time horizon; tau = [-1,+1]
p.D   = LGL_Dmatrix(p.tau); % For defect constraints
p.w   = LGL_weights(p.tau); % For gaussian quadrature

%% Initial Conditions
X0 = zeros(p.nt*(p.ns+p.nu),1); % Initial Guess (all zeros)

%% Solution

options = optimoptions(@fmincon,'Display','final','MaxFunEvals',1e5);
[sol,val,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x) objfun(x,p),X0,[],[],[],[],[],[],@(x) constraints(x,p),options);

%% Results

t = (p.tau*(p.tf-p.t0)+(p.tf+p.t0))/2;
x = sol(p.xi);
z = sol(p.zi);
u = sol(p.ui);

lamx  = lambda.eqnonlin(p.xi)./p.w;
lamz  = lambda.eqnonlin(p.zi)./p.w;

figure(1);hold on;
subplot(221);
plot(x,z,'r.-');grid
xlabel('x');ylabel('z');
subplot(222);
plot(t,lamx,'r.-',t,lamz,'b.-');grid
legend('lamx','lamz','Location','Best');
xlabel('t');ylabel('costates');
subplot(223);
plot(t,x,'r.-',t,z,'b.-');grid
legend('x','z','Location','Best');
xlabel('t');ylabel('states');
subplot(224);
plot(t,tan(u),'r.-',t,lamz./lamx,'b--');grid
legend('tan(\gamma)','lamz/lamx','Location','Best');
xlabel('t');ylabel('control');

%%
figure(2)
title('Costates');
subplot(211);hold on;grid on
plot(t,-lamx,'r.-')
plot(t,-1*ones(length(t),1),'bo');
ylabel('lamx');
legend('Numerical','Analytical');
subplot(212);hold on;grid on
plot(t,-lamz,'r.-')
plot(t,-p.K+p.K.*t,'bo');
ylabel('lamz');
legend('Numerical','Analytical');

%% Objective Function 

function f = objfun(x,p)
f  = -x(p.xi(end));
end

%% Constraints

function [c,ceq] = constraints(x,p)

xx = x(p.xi);
zz = x(p.zi);
u  = x(p.ui);

%% Inequality Constraints

c = [u-pi;-pi-u]; % -pi <= u <= pi

%% Equality Constraints

Y = [xx,zz]; % States
F = (p.tf-p.t0)/2*[p.V*cos(u)+p.K.*zz , p.V*sin(u)]; % RHS of ODE

ceq1 = p.D*Y-F;    % Defect Constraints
ceq2 = xx(1)-p.x0; % Initial Value Constraint on x 
ceq3 = zz(1)-p.z0; % Initial Value Constraint on z

ceq = [ceq1(:);ceq2;ceq3;];

end