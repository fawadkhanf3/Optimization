%% Brachistochrone Problem (Solved By LGL Pseudospectral Method)

clear;close all;clc
format long g
set(0,'DefaultLineLineWidth',2);

%% Problem Parameters

p.ns = 2; % Number of States
p.nu = 1; % Number of Controls

p.t0 = 0; % Initial time

p.x0 = 0; % Initial condition for x
p.y0 = 0; % Initial condition for y
p.xf = 1; % Final condition for x

p.g  = 9.81; 

p.nt = 10; % Number of Node Points

% Discretized variable indices in X = [x,z,u];
p.xi = 1:p.nt;
p.yi = p.nt+1:2*p.nt;
p.ui = 2*p.nt+1:3*p.nt;
p.tfi = 3*p.nt+1;

%% LGL Nodes, Differential Matrix and Weights

p.tau = LGL_nodes(p.nt-1);  % Scaled time horizon; tau = [-1,+1]
p.D   = LGL_Dmatrix(p.tau); % For defect constraints
p.w   = LGL_weights(p.tau); % For gaussian quadrature

%% Initial Conditions
% X0 = zeros(p.nt*(p.ns+p.nu),1); % Initial Guess (all zeros)
X0(p.xi)  = p.x0;
X0(p.yi)  = p.y0;
X0(p.ui)  = 1.5;
X0(p.tfi) = 0.5;

X0 = X0';

%% Solution

options = optimoptions(@fmincon,'Display','iter',...
    'MaxFunEvals',1e5,'StepTolerance',1e-20,'MaxIterations',5000);
[sol,val,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x) objfun(x,p),X0,[],[],[],[],[],[],@(x) constraints(x,p),options);

%% Results

tf = sol(p.tfi);
t = (p.tau*(tf-p.t0)+(tf+p.t0))/2;
x = sol(p.xi);
y = sol(p.yi);
u = sol(p.ui);

lamx  = -lambda.eqnonlin(p.xi)./p.w;
lamy  = -lambda.eqnonlin(p.yi)./p.w;

figure(1);hold on;grid on
plot(x,y,'r.-');
xlabel('x');ylabel('y');
title('States');

figure(2);hold on;grid on
plot(t,u,'r.-');
xlabel('t');ylabel('\theta');
title('Control');

figure(3);hold on;grid on
plot(t,lamx,'r.-');
plot(t,-1/2*sqrt(pi/p.g)*ones(length(t),1),'b--');
xlabel('t');ylabel('lamx')
legend('Numerical','Analytical','Location','Best');

figure(4);hold on;grid on
plot(t,lamy,'r.-');
plot(t,lamx.*cot((sqrt(pi*p.g)/2).*t),'b--');
xlabel('t');ylabel('lamy')
legend('Numerical','Analytical','Location','Best');

%% Objective Function 

function f = objfun(x,p)
f  = x(p.tfi);
end

%% Constraints

function [c,ceq] = constraints(x,p)

xx = x(p.xi);
yy = x(p.yi);
u  = x(p.ui);
tf = x(p.tfi);

%% Inequality Constraints

c = -tf;%[u-pi;-pi-u]; % -pi <= u <= pi

%% Equality Constraints

Y = [xx,yy]; % States
F = (tf-p.t0)/2*[real(sqrt(2*p.g.*yy)).*cos(u),real(sqrt(2*p.g.*yy)).*sin(u)]; % RHS of ODE

ceq1 = p.D*Y-F;    % Defect Constraints
ceq2 = xx(1)-p.x0; % Initial Value Constraint on x 
ceq3 = yy(1)-p.y0; % Initial Value Constraint on y
ceq4 = xx(end)-p.xf; % Final Value Constraint on x
ceq = [ceq1(:);ceq2;ceq3;ceq4];

end