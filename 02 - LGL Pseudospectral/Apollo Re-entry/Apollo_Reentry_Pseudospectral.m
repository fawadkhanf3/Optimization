%% Apollo Reentry Problem (Solved By LGL Pseudospectral Method)

% Example 2 from Direct & Indirect Methods for Trajectory Optimization
% by Bulirsh and Stryk

clear;close all;clc
format long g
set(0,'DefaultLineLineWidth',2);

%% Problem Parameters

p.ns = 5; % Number of States
p.nu = 1; % Number of Controls

p.t0 = 0; % Initial time

%%

p.R    = 209.0352;
p.g0   = 3.2172e-4;
p.Sbym = 50000;
p.beta = 1/0.235;
p.c    = 20;
p.N    = 4;

p.rho0 = 2.3769e-3;
p.rho  = @(zi) p.rho0.*exp(-p.beta.*p.R.*zi);

p.cL0  = -0.505;
p.cD0  = 0.88;
p.cDL  = 0.52;

p.cD   = @(u) p.cD0 + p.cDL.*cos(u);
p.cL   = @(u) p.cL0.*sin(u);

p.v0     = 0.35;
p.zi0    = 4/p.R;
p.zeta0  = 0;
p.gamma0 = -5.75*pi/180;
p.q0     = 0;

p.vf     = 0.0165;
p.zif    = 0.75530/p.R;
p.zetaf  = 51.6912;

p.nt = 64; % Number of Node Points

% Discretized variable indices in X = [x,z,u];
p.vi     = 1:p.nt;
p.zii    = p.nt+1:2*p.nt;
p.zetai  = 2*p.nt+1:3*p.nt;
p.gammai = 3*p.nt+1:4*p.nt;
p.qi     = 4*p.nt+1:5*p.nt;
p.ui     = 5*p.nt+1:6*p.nt;
p.tfi    = 6*p.nt+1;

%% LGL Nodes, Differential Matrix and Weights

p.tau = LGL_nodes(p.nt-1);  % Scaled time horizon; tau = [-1,+1]
p.D   = LGL_Dmatrix(p.tau); % For defect constraints
p.w   = LGL_weights(p.tau); % For gaussian quadrature

%% Initial Conditions
% X0 = zeros(p.nt*(p.ns+p.nu),1); % Initial Guess (all zeros)

X0(p.vi)     = p.v0;
X0(p.zii)    = p.zi0;
X0(p.zetai)  = p.zeta0;
X0(p.gammai) = p.gamma0;
X0(p.qi)     = p.q0;
X0(p.ui)     = 90*pi/180;
X0(p.tfi)    = 400;

X0 = X0';

%% Solution

options = optimoptions(@fmincon,'Display','iter',...
    'MaxFunEvals',1e6,'StepTolerance',1e-20,'MaxIterations',5000);
[sol,val,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x) objfun(x,p),X0,[],[],[],[],[],[],@(x) constraints(x,p),options);

%% Results

tf = sol(p.tfi);

t = (p.tau*(tf-p.t0)+(tf+p.t0))/2;

v     = sol(p.vi);
zi    = sol(p.zii);
zeta  = sol(p.zetai);
gamma = sol(p.gammai);
q     = sol(p.qi);
u     = sol(p.ui);

lamv     = -lambda.eqnonlin(p.vi)./p.w;
lamzi    = -lambda.eqnonlin(p.zii)./p.w;
lamzeta  = -lambda.eqnonlin(p.zetai)./p.w;
lamgamma = -lambda.eqnonlin(p.gammai)./p.w;
lamq     = -lambda.eqnonlin(p.qi)./p.w;

figure(1);hold on;grid on;
plot(t,zi,'r.-');grid;shg

%% Plots

%% Objective Function 

function f = objfun(x,p)
q = x(p.qi);
f = q(end);
end

%% Constraints

function [c,ceq] = constraints(x,p)

v     = x(p.vi);
zi    = x(p.zii);
zeta  = x(p.zetai);
gamma = x(p.gammai);
q     = x(p.qi);
u     = x(p.ui);
tf    = x(p.tfi);

%% Inequality Constraints

c = [u-pi;-pi-u]; % -pi <= E <= pi

%% Equality Constraints

Y = [v,zi,zeta,gamma,q]; % States

dv     = -(1/2*p.Sbym).*(p.rho(zi)).*(v.^2).*(p.cD(u))- (p.g0.*sin(gamma))./((1+zi).^2);

dzi    = (v./p.R) .* sin(gamma);

dzeta  = (v./(1+zi)).*cos(gamma);

dgamma = (1/2*p.Sbym).*(p.rho(zi)).*(v).*(p.cL(u)) + (v.*cos(gamma))./(p.R.*(1+zi)) - ...
         (p.g0.*cos(gamma))./(v.*(1+zi).^2);

dq = (p.c).*(v.^3).*(sqrt((p.rho(zi))./p.N));

F = (tf-p.t0)/2*[dv,dzi,dzeta,dgamma,dq]; % RHS of ODE

ceq1 = p.D*Y-F;    % Defect Constraints

ceq2 = v(1)-p.v0;
ceq3 = zi(1)-p.zi0;
ceq4 = zeta(1)-p.zeta0;
ceq5 = gamma(1)-p.gamma0;
ceq6 = q(1)-p.q0;
ceq7 = v(end)-p.vf;
ceq8 = zi(end)-p.zif;
ceq9 = zeta(end)-p.zetaf;

ceq = [ceq1(:);ceq2;ceq3;ceq4;ceq5;ceq6;ceq7;ceq8;ceq9];

end