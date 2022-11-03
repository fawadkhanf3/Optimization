%% Orbit Transfer Problem (Solved By LGL Pseudospectral Method)

% Example 1 from Costate Estimation by a Legendre Pseudospectral Method 
% Fariba Fahroo & Micheal Ross

clear;close all;clc
format long g
set(0,'DefaultLineLineWidth',2);

%% Problem Parameters

p.ns = 4; % Number of States
p.nu = 1; % Number of Controls

p.t0 = 0; % Initial time
p.tf = 50; % Final time

p.r0     = 1.1;
p.theta0 = 0;
p.u0     = 0;
p.v0     = 1/sqrt(1.1);

p.A  = 0.01;

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

%% Initial Conditions
% X0 = zeros(p.nt*(p.ns+p.nu),1); % Initial Guess (all zeros)

X0(p.ri)     = p.r0;
X0(p.thetai) = p.theta0;
X0(p.ui)     = p.u0;
X0(p.vi)     = p.v0;
X0(p.Ei)     = 0.05;

X0 = X0';

%% Solution

options = optimoptions(@fmincon,'Display','iter',...
    'MaxFunEvals',1e6,'StepTolerance',1e-20,'MaxIterations',5000);
[sol,val,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x) objfun(x,p),X0,[],[],[],[],[],[],@(x) constraints(x,p),options);

save Solution.mat sol p
%% Results

t = (p.tau*(p.tf-p.t0)+(p.tf+p.t0))/2;

r     = sol(p.ri);
theta = sol(p.thetai);
u     = sol(p.ui);
v     = sol(p.vi);
E     = sol(p.Ei);

lamr     = -lambda.eqnonlin(p.ri)./p.w;
lamtheta = -lambda.eqnonlin(p.thetai)./p.w;
lamu     = -lambda.eqnonlin(p.ui)./p.w;
lamv     = -lambda.eqnonlin(p.vi)./p.w;

%% Plots
figure(1);hold on;grid on
plot(r.*cos(theta),r.*sin(theta),'ro-');

figure(2);hold on;grid on
plot(t,u,'ro-',t,v,'bo-');
xlabel('t [sec]');
ylabel('Optimal Velocities');
legend('u','v');

figure(3);hold on;grid on
plot(t,lamr,'ro-',t,lamu,'bo-',t,lamv,'ko-');
xlabel('t [sec]');
ylabel('Costates');
legend('lamr','lamu','lamv');

figure(4);hold on;grid on
plot(t,tan(E),'r.-',t,lamu./lamv,'bo');
xlabel('t [sec]');
ylabel('Control Angle');
legend('Direct','Costate Estimates');

nInter = 2000;
tArray = linspace(p.t0,p.tf,nInter);
r2     = LagrangeInter(t',r',tArray);
theta2 = LagrangeInter(t',theta',tArray);
figure(5);hold on;grid on;box on;
plot(r2.*cos(theta2),r2.*sin(theta2),'r.-');shg

%%

r3     = interp1(t,r,tArray);
theta3 = interp1(t,theta,tArray);
figure(6);hold on;grid on;box on;
plot(r3.*cos(theta3),r3.*sin(theta3),'r.-');shg

%% Objective Function 

function f = objfun(x,p)

r = x(p.ri);
u = x(p.ui);
v = x(p.vi);

f  = -(1/2*(u(end)^2 + v(end)^2)-1/r(end));
end

%% Constraints

function [c,ceq] = constraints(x,p)

r     = x(p.ri);
theta = x(p.thetai);
u     = x(p.ui);
v     = x(p.vi);
E     = x(p.Ei);

%% Inequality Constraints

c = [E-pi;-pi-E]; % -pi <= E <= pi

%% Equality Constraints

Y = [r,theta,u,v]; % States

dr     = u;
dtheta = v./r;
du     = (v.^2)./r - 1./(r.^2) + p.A.*sin(E);
dv     = -(u.*v)./r + p.A.*cos(E);

F = (p.tf-p.t0)/2*[dr,dtheta,du,dv]; % RHS of ODE

ceq1 = p.D*Y-F;    % Defect Constraints

ceq2 = r(1)-p.r0;
ceq3 = theta(1)-p.theta0;
ceq4 = u(1)-p.u0;
ceq5 = v(1)-p.v0;

ceq = [ceq1(:);ceq2;ceq3;ceq4;ceq5];

end