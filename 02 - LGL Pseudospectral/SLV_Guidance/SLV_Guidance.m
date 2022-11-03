%% SLV Guidance (Solved By LGL Pseudospectral Method)

clear;close all;clc
format long g
set(0,'DefaultLineLineWidth',2);

%% Problem Parameters

p.ns = 4; % Number of States
p.nu = 1; % Number of Controls

p.t0 = 0; % Initial time

p.G  = 6.67e-11;
p.M  = 5.972e24;
p.RE = 6371e3;
p.g0 = 9.81;

p.Thr  = 14.709975e3;
p.m0   = 3630;
p.mdot = 4.765;

p.x0  = 0; 
p.y0  = 7085e03; 
p.vx0 = 5575;
p.vy0 = 843.264;

p.rf = 824e3+p.RE;
p.Vf = sqrt(p.G*p.M/p.rf);
p.gammaf = 0.0;

p.nt = 40; % Number of Node Points

% Discretized variable indices in X = [x,y,vx,vy,u,tf];
p.xi  = 1:p.nt;
p.yi  = p.nt+1:2*p.nt;
p.vxi = 2*p.nt+1:3*p.nt;
p.vyi = 3*p.nt+1:4*p.nt;
p.ui  = 4*p.nt+1:5*p.nt;
p.tfi = 5*p.nt+1;

%% LGL Nodes, Differential Matrix and Weights

p.tau = LGL_nodes(p.nt-1);  % Scaled time horizon; tau = [-1,+1]
p.D   = LGL_Dmatrix(p.tau); % For defect constraints
p.w   = LGL_weights(p.tau); % For gaussian quadrature

%% Initial Conditions

% X0(p.xi)  = p.x0;
% X0(p.yi)  = p.y0;
% X0(p.vxi) = p.vx0;
% X0(p.vyi) = p.vy0;
% X0(p.ui)  = -10*pi/180;
% X0(p.tfi) = 350;
% 
% X0 = X0';

Guess = load('Guess.mat');
tg   = (p.tau*(Guess.t(end)-p.t0)+(Guess.t(end)+p.t0))/2;

X0(p.xi,1)  = interp1(Guess.t,Guess.x,tg);
X0(p.yi,1)  = interp1(Guess.t,Guess.y,tg);
X0(p.vxi,1) = interp1(Guess.t,Guess.vx,tg);
X0(p.vyi,1) = interp1(Guess.t,Guess.vy,tg);
X0(p.ui,1)  = interp1(Guess.t,Guess.u,tg);
X0(p.tfi,1) = Guess.t(end);


%% Solution

options = optimoptions(@fmincon,'Display','iter',...
    'MaxFunEvals',1e6,'MaxIterations',5000,'StepTolerance',1e-20);

[sol,val,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x) objfun(x,p),X0,[],[],[],[],[],[],@(x) constraints(x,p),options);

%% Results

tf = sol(p.tfi);

t = (p.tau*(tf-p.t0)+(tf+p.t0))/2;

x  = sol(p.xi);
y  = sol(p.yi);
vx = sol(p.vxi);
vy = sol(p.vyi);
u  = sol(p.ui);

save Converged.mat sol lambda p tf t x y vx vy u options X0

% plt2

%% Objective Function 

function f = objfun(x,p)
f  = x(p.tfi);
end

%% Constraints

function [c,ceq] = constraints(x,p)

xx = x(p.xi);
yy = x(p.yi);
vx = x(p.vxi);
vy = x(p.vyi);
u  = x(p.ui);
tf = x(p.tfi);

%% Inequality Constraints

c = [u-pi/2 -pi/2-u]; % -pi <= u <= pi

%% Equality Constraints

Y = [xx,yy,vx,vy]; % States

t = (p.tau*(tf-p.t0)+(tf+p.t0))/2;
r = sqrt(xx.^2+yy.^2);

mass = p.m0-p.mdot.*t;
aT   = p.Thr./mass;
dvx  = aT.*cos(u) - (p.G*p.M.*xx)./(r.^3);
dvy  = aT.*sin(u) - (p.G*p.M.*yy)./(r.^3);

F = (tf-p.t0)/2*[vx,vy,dvx,dvy]; % RHS of ODE

ceq1 = p.D*Y-F;    % Defect Constraints

% Initial Constraints
ceq2 = xx(1)-p.x0;
ceq3 = yy(1)-p.y0; 
ceq4 = vx(1)-p.vx0;
ceq5 = vy(1)-p.vy0;

% Terminal Constraints
rfinal = [xx(end);yy(end)];
vfinal = [vx(end);vy(end)];

ceq6 = norm(rfinal)-p.rf;
ceq7 = norm(vfinal)-p.Vf;
ceq8 = dot(rfinal,vfinal)-p.gammaf;

ceq = [ceq1(:);ceq2;ceq3;ceq4;ceq5;ceq6;ceq7;ceq8];

end