%% Nonlinear Simulation of Inverted Pendulum on a Cart

clear; close all; clc;
format short g

s = tf('s');

%% Inverted Pendulum
M = 1.0; 
m = 0.25;
l = 1.0;
b = 0.05;
g = 9.81;

%% State Space
A = [0 0 1 0;....
     0 0 0 1;...
     0 -(m*g/M) -b/M 0;...
     0 (2*(M + m)*g)/(M*l) 2*b/(M*l) 0];
B = [0; 0; 1/M; -2/(M*l)];
C = [1 0 0 0;0 1 0 0];
D = [0;0];
sys = ss(A,B,C,D);

%% Controller Design

Transfer_function = tf(sys);

theta = Transfer_function(2,:);

options = optimoptions('particleswarm',...
    'SwarmSize',20,'HybridFcn',@fmincon);
rng(8)  % For reproducibility
nvars = 7;
bounds = [1 2;  % [Mmin, Mmax] <=> Wp
          0 1;  % [Amin, Amax] <=> Wp 
          0 30; % [wbmin, wbmax] <=> Wp
          1 2;  % [Mmin, Mmax] <=> Wt
          0 1;  % [Amin, Amax] <=> Wt
          0 30; % [wb_min, wb_max] <=> Wt
          0 5]; % [cmin, cmax] <=> Wu
lb = bounds(:,1);
ub = bounds(:,2);
    
x = particleswarm(@(x) obj_fun(x,theta),nvars,lb,ub,options);

%%
wp.M  = x(1);
wp.A  = x(2);
wp.wB = x(3);

Wp = tf([1/wp.M wp.wB],[1 wp.wB*wp.A]);

wt.M  = x(4);
wt.A  = x(5);
wt.wB = x(6);

Wt = tf([1/wt.M wt.wB],[1 wt.wB*wt.A]);

Wu = tf(x(7));

[K,CL,GAM,INFO] = mixsyn(theta,Wp,Wu,Wt,'Display','On');

L = K*theta;
S = 1/(1+L);
T = 1-S;

% set(0,'DefaultLineLineWidth',2);
% figure(1);grid on;box on;
% sigma(S,'b',1/Wp,'r--',{.1,3000});
% legend('S','1/Wp','Location','Best');
% figure(2);grid on;box on;
% sigma(T,'b',1/Wt,'r--',{.1,3000});
% legend('T','1/Wt','Location','Best');
% return

%% Nominal Analysis

sim('NonlinearSim.slx')

figure(1),clf
plot(time, angle,'b','linewidth',2),hold on
xlabel('time (seconds)');ylabel('angle (degrees)');grid on;
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times');

figure(2),clf
plot(time, position,'b','linewidth',2),hold on
xlabel('time (seconds)');ylabel('position (meters)');grid on;
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times');

%%

aug_sys = feedback(K*sys,1,1,2);
aug_sys_tf = tf(aug_sys);
position = aug_sys_tf(1,:);

options = optimoptions('particleswarm','SwarmSize',20,'HybridFcn',@fmincon);
rng(8)  % For reproducibility
nvars = 7;
bounds = [1 2;  % [Mmin, Mmax] <=> Wp
          0 1;  % [Amin, Amax] <=> Wp 
          0 30; % [wbmin, wbmax] <=> Wp
          1 2;  % [Mmin, Mmax] <=> Wt
          0 1;  % [Amin, Amax] <=> Wt
          0 30; % [wb_min, wb_max] <=> Wt
          0 5]; % [cmin, cmax] <=> Wu
lb = bounds(:,1);
ub = bounds(:,2);
    
x = particleswarm(@(x) obj_fun(x,position),nvars,lb,ub,options);

%%
wp.M  = x(1);
wp.A  = x(2);
wp.wB = x(3);

Wp = tf([1/wp.M wp.wB],[1 wp.wB*wp.A]);

wt.M  = x(4);
wt.A  = x(5);
wt.wB = x(6);

Wt = tf([1/wt.M wt.wB],[1 wt.wB*wt.A]);

Wu = tf(x(7));

[K2,CL,GAM,INFO] = mixsyn(position,Wp,Wu,Wt,'Display','On');

%% Nominal Analysis

sim('NonlinearSim2.slx')

figure(1)
plot(time, angle,'r','linewidth',2)

figure(2)
plot(time, position,'r','linewidth',2)
