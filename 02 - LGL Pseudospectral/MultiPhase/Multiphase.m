%% Multiphase Problem (Solved By LGL Pseudospectral Method)

clear;close all;clc
warning off
format long g
set(0,'DefaultLineLineWidth',2);

%% Problem Parameters

p.ns = 5; % Number of States
p.nu = 1; % Number of Controls

p.x0  = 36.022e3;
p.y0  = 60.708e3;
p.vx0 = 1.052e3;
p.vy0 = 1.060e3;
p.m0  = 76501;

p.xf  = 0.0e3;
p.yf  = 0.0e3;
p.vxf = 0.0e3;
p.vyf = 0.5e3;

p.TLO  = 5886e3;
p.Isp  = 282;
p.d    = 3.66;
p.Cd   = 0.75;
p.h0   = 7500;
p.g0   = 9.80665;
p.rho0 = 1.225;

p.t0   = 0;
p.t1   = 40.8;
p.mdry = 25600;

p.T1   = 1/3*p.TLO;
p.T2   = 0;
p.T3   = 1/9*p.TLO;

p.k1   = 1;
p.k2   = 0;
p.k3   = 1;

p.nt = 16; % Number of Node Points

% Discretized variable indices in X = [x,y,vx,vy,u,tf];
p.x1i  = 1:p.nt;
p.y1i  = p.nt+1:2*p.nt;
p.vx1i = 2*p.nt+1:3*p.nt;
p.vy1i = 3*p.nt+1:4*p.nt;
p.m1i  = 4*p.nt+1:5*p.nt;
p.u1i  = 5*p.nt+1:6*p.nt;

p.x2i  = 6*p.nt+1:7*p.nt;
p.y2i  = 7*p.nt+1:8*p.nt;
p.vx2i = 8*p.nt+1:9*p.nt;
p.vy2i = 9*p.nt+1:10*p.nt;
p.m2i  = 10*p.nt+1:11*p.nt;
p.u2i  = 11*p.nt+1:12*p.nt;

p.x3i  = 12*p.nt+1:13*p.nt;
p.y3i  = 14*p.nt+1:15*p.nt;
p.vx3i = 15*p.nt+1:16*p.nt;
p.vy3i = 16*p.nt+1:17*p.nt;
p.m3i  = 17*p.nt+1:18*p.nt;
p.u3i  = 18*p.nt+1:19*p.nt;

p.tf2i = 19*p.nt+1; 
p.tf3i = 19*p.nt+2;

%% LGL Nodes, Differential Matrix and Weights

p.tau = LGL_nodes(p.nt-1);  % Scaled time horizon; tau = [-1,+1]
p.D   = LGL_Dmatrix(p.tau); % For defect constraints
p.w   = LGL_weights(p.tau); % For gaussian quadrature

%% Initial Conditions

X0(p.x1i)  = p.x0;
X0(p.y1i)  = p.y0;
X0(p.vx1i) = p.vx0;
X0(p.vy1i) = p.vy0;
X0(p.m1i)  = p.m0;
X0(p.u1i)  = -10*pi/180;

X0(p.x2i)  = p.x0;
X0(p.y2i)  = p.y0;
X0(p.vx2i) = p.vx0;
X0(p.vy2i) = p.vy0;
X0(p.m2i)  = p.m0-p.mdry;
X0(p.u2i)  = -10*pi/180;

X0(p.x3i)  = p.x0;
X0(p.y3i)  = p.y0;
X0(p.vx3i) = p.vx0;
X0(p.vy3i) = p.vy0;
X0(p.m3i)  = p.m0-p.mdry;
X0(p.u3i)  = -10*pi/180;

X0(p.tf2i) = 220;
X0(p.tf3i) = 305;

X0 = X0';

%% Solution

options = optimoptions(@fmincon,'Display','iter',...
    'MaxFunEvals',1e6,'MaxIterations',5000,'StepTolerance',1e-20);

[sol,val,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x) objfun(x,p),X0,[],[],[],[],[],[],@(x) constraints(x,p),options);

%% Results

% tf = sol(p.tfi);
% 
% t = (p.tau*(tf-p.t0)+(tf+p.t0))/2;
% 
% x  = sol(p.xi);
% y  = sol(p.yi);
% vx = sol(p.vxi);
% vy = sol(p.vyi);
% u  = sol(p.ui);
% 
% save Converged.mat sol lambda p tf t x y vx vy u options X0

%% Objective Function 

function f = objfun(x,p)

mass = x(p.m3i);
f    = -mass(end);

end

%% Constraints

function [c,ceq] = constraints(x,p)

x1  = x(p.x1i);
y1  = x(p.y1i);
vx1 = x(p.vx1i);
vy1 = x(p.vy1i);
m1  = x(p.m1i);
u1  = x(p.u1i);

x2  = x(p.x2i);
y2  = x(p.y2i);
vx2 = x(p.vx2i);
vy2 = x(p.vy2i);
m2  = x(p.m2i);
u2  = x(p.u2i);

x3  = x(p.x3i);
y3  = x(p.y3i);
vx3 = x(p.vx3i);
vy3 = x(p.vy3i);
m3  = x(p.m3i);
u3  = x(p.u3i);

tf1 = p.t1;
tf2 = x(p.tf2i);
tf3 = x(p.tf3i);

%% Equality Constraints

t1 = (p.tau*(tf1-p.t0)+(tf1+p.t0))/2;
t2 = (p.tau*(tf2-tf1)+(tf2+tf1))/2;
t3 = (p.tau*(tf3-tf2)+(tf3+tf2))/2;

[dx1,dy1,dvx1,dvy1,dm1] = derivs(t1,x1,y1,vx1,vy1,m1,u1,p.T1,p.k1,p);
[dx2,dy2,dvx2,dvy2,dm2] = derivs(t2,x2,y2,vx2,vy2,m2,u2,p.T2,p.k2,p);
[dx3,dy3,dvx3,dvy3,dm3] = derivs(t3,x3,y3,vx3,vy3,m3,u3,p.T3,p.k3,p);

Y1 = [x1,y1,vx1,vy1,m1];
Y2 = [x2,y2,vx2,vy2,m2];
Y3 = [x3,y3,vx3,vy3,m3]; % States

F1 = (tf1-p.t0)/2*[dx1,dy1,dvx1,dvy1,dm1];
F2 = (tf2-tf1)/2*[dx2,dy2,dvx2,dvy2,dm2];
F3 = (tf3-tf2)/2*[dx3,dy3,dvx3,dvy3,dm3];

ceq1 = p.D*Y1-F1;    % Defect Constraints
ceq2 = p.D*Y2-F2;    % Defect Constraints
ceq3 = p.D*Y3-F3;    % Defect Constraints

ceq4 = x1(1)  - p.x0;
ceq5 = y1(1)  - p.y0;
ceq6 = vx1(1) - p.vx0;
ceq7 = vy1(1) - p.vy0;
ceq8 = m1(1)  - p.m0;

ceq9  = x3(end) - p.xf;
ceq10 = y3(end) - p.yf;
ceq11 = vx3(end) - p.vxf;

ceq12 = x2(1)  - x1(end);
ceq13 = y2(1)  - y1(end);
ceq14 = vx2(1) - vx1(end);
ceq15 = vy2(1) - vy1(end);
ceq16 = m2(1)  - m1(end) - p.mdry;

ceq17 = x3(1)  - x2(end);
ceq18 = y3(1)  - y2(end);
ceq19 = vx3(1) - vx2(end);
ceq20 = vy3(1) - vy2(end);
ceq21 = m3(1) - m2(end);

c = [-p.vyf-vy3(end); vy3(end)-p.vyf];

ceq = [ceq1(:);ceq2(:);ceq3(:);ceq4;ceq5;ceq6;ceq7;ceq8;...
    ceq9;ceq10;ceq11;ceq12;ceq13;ceq14;ceq15;ceq16;ceq17;...
    ceq18;ceq19;ceq20;ceq21];



end

function [dx,dy,dvx,dvy,dm] = derivs(t,x,y,vx,vy,m,u,T,k,p)

dx = vx;
dy = vy;

v   = sqrt(vx.^2+vy.^2);

cgam = vx./v;
sgam = vy./v;

D   = drag(y,v,p);

dvx = T./m .* k .* cos(u) - D./m .* cgam;
dvy = T./m .* k .* sin(u) - D./m .* sgam - p.g0;

dm = -T./p.Isp/p.g0 * k * ones(size(m));

end

function D = drag(y,v,p)

rho0 = p.rho0;
h0   = p.h0;
Cd   = p.Cd;
d    = p.d;
A    = pi*d^2/4;
v2   = v.^2;

term = exp(-y./h0);
D = 0.5*rho0.*term.*Cd.*A.*v2;

end