%% Optimal Control Problem (Trajectory Optimization)

%------------------------------------------------------------%
function solution = Optimal_Control_Problem()

format long

% Universal Parameters
global G M RE g
G = 6.67e-11; % Gravitational Constant 
M = 5.972e24; % Mass of Earth
RE = 6371e03; % Radius of Earth (m)
g = 9.81; % Gravitational Acceleration (m/s^2)

% Mission Specifications
global h rf gamma Vc
h = 824e03; % Orbit Altitude (m)
rf = h + RE; % Orbit Radius (m)
gamma = 0; % Flight Path Angle (degrees)
Vc = 7443.3; % Orbital Velocity (m/s)

% SLV Specifications
global T mo mdot 
T = 14.709975e03; % Thrust (N)
mo = 3630; % Propellant Mass (kg)
mdot = 4.765; % Propellant Consumption Rate(kg/sec)

% Initial Conditions (At the Last Stage Ignition)
global xo yo vxo vyo
xo = 0; % Initial x position (m)
yo = 7085e03 ; % Initial y position (m)
vxo = 5575; % Initial x velocity (m/s)
vyo = 843.264; % Initial y velocity (m/s)

% Final Conditions
global Yf Vxf Vyf
Yf = rf; % Final Altitude
Vxf = Vc; % Final Downrange Velocity
Vyf = 0; % Final Vertical Velocity

% Initial Guesses
global tau Nt
to = 0; % Initial Time
Nt = 120; % Number of divisions

tau = linspace(0,1,Nt); % Non-dimensional Time Vector

yinit = [ones(1,Nt)*xo;ones(1,Nt)*yo;ones(1,Nt)*vxo;ones(1,Nt)*vyo;... % Initial Guess for States
         ones(1,Nt)*4e-4;-ones(1,Nt)*6e-4;-ones(1,Nt)*1;ones(1,Nt)*0.2]; % Initial Guess for Costates
        % [ x_guess; y_guess; vx_guess; vy_guess;
        %   lambda_x_guess; lambda_y_guess; lambda_vx_guess; lambda_vy_guess]

tf_guess = 700; % Guess for Final Time

% Solution Initialization
solinit.x = tau;
solinit.y = yinit;
solinit.parameters = tf_guess;
   % solinit = bvpinit(tau,yinit,tf_guess);

tol = 1e-10; % Tolerance

% Solver Options
options = bvpset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol tol tol],'Nmax',2000);
    % RelTol = Relative Tolerance
    % AbsTol = Absolute Tolerance
    % Nmax = Maximum Number of Mesh Points

% Solution
sol = bvp4c(@orbit_odes,@orbit_bcs,solinit,options);
    % @orbit_odes = Calls the function with Differential Equations
    % @orbit_bcs = Calls the funtion with the Boundary Conditions

% Final Time Value
tf = sol.parameters;
disp('Optimal Burn Time for the Last Stage = ')
disp(tf)

% Residual Values Obtained to achieve Final Conditions
[res] = orbit_bcs(sol.y(:,1),sol.y(:,end),sol.parameters);

% Solution at all times in time vector tau
Z = deval(sol,tau);

% Dimensional Time
time = sol.x*tf;

% Individual Solutions
x_sol = sol.y(1,:);
y_sol = sol.y(2,:);
vx_sol = sol.y(3,:);
vy_sol = sol.y(4,:);
lambda_x_sol = sol.y(5,:);
lambda_y_sol = sol.y(6,:);
lambda_vx_sol = sol.y(7,:);
lambda_vy_sol = sol.y(8,:);
radius = sqrt((sol.y(1,:).^2)+(sol.y(2,:).^2));
velocity = sqrt((sol.y(3,:).^2)+(sol.y(4,:).^2));
for i=1:length(sol.x)
    xx1 = [sol.y(1,i);sol.y(2,i)];
    vv1 = [sol.y(3,i);sol.y(4,i)];
    gama(i) = acos(dot(xx1,vv1)/(norm(xx1)*norm(vv1)));
end
control_angle = atan2(-sol.y(8,:),-sol.y(7,:));

solution = {time,x_sol,y_sol,vx_sol,vy_sol,control_angle};

% Differential Equations for States & Co-states
function dydx = orbit_odes(tau,x1,tf)

% State Variables
x = x1(1); % Horizontal Position
y = x1(2); % Vertical Position
vx = x1(3); % Horizontal Velocity
vy = x1(4); % Vertical Velocity

% Co-state Variables
lambda_x = x1(5);
lambda_y = x1(6);
lambda_vx = x1(7);
lambda_vy = x1(8);

global G M mo mdot T m

m = (mo - mdot*tau*tf); % Mass w.r.t time
aT = T/m; % Acceleration (variable)

r = sqrt((x^2)+(y^2)); % Radius
lambda_v = sqrt((lambda_vx^2)+(lambda_vy^2));
xx = [x;y];
lamv = [lambda_vx;lambda_vy];
ux = -lambda_vx/norm(lamv);
uy = -lambda_vy/norm(lamv);
constant = 3*dot(xx,lamv)/(r^2);

% Differential Equations for States
d_x = vx;
d_y = vy;
d_vx =  aT*ux - (G*M*x/(r^3));
d_vy =  aT*uy - (G*M*y/(r^3));

% Differential Equations for Co-states
d_lambda_x = (G*M/(r^3))*(lambda_vx - constant*x);
d_lambda_y = (G*M/(r^3))*(lambda_vy - constant*y);
d_lambda_vx = -lambda_x;
d_lambda_vy = -lambda_y;

dydx = tf*[d_x;d_y;d_vx;d_vy;d_lambda_x;d_lambda_y;d_lambda_vx;d_lambda_vy];
return;
%------------------------------------------------------------%

% Boundary Conditions (Initial & Final Conditions)

function res = orbit_bcs(x1,xb,tf)

% Initial States
x = x1(1); 
y = x1(2);
vx = x1(3);
vy = x1(4);
lambda_x = x1(5);
lambda_y = x1(6);
lambda_vx = x1(7);
lambda_vy = x1(8);

% Final States
xf = xb(1); 
yf = xb(2);
vxf = xb(3);
vyf = xb(4);
lambda_xf = xb(5);
lambda_yf = xb(6);
lambda_vxf = xb(7);
lambda_vyf = xb(8);

global G M rf xo yo vxo vyo
rfinal = [xf;yf];
vfinal = [vxf;vyf];

% Boundary Conditions

% Initial Conditions
b1 = x - xo;
b2 = y - yo;
b3 = vx - vxo;
b4 = vy - vyo;

% Optimality Condition
b5 = sqrt(lambda_vx*lambda_vx+lambda_vy*lambda_vy)-1;

% Final End Conditions
b6 = norm(rfinal)-rf;%sqrt((xf^2)+(yf^2))-rf;
b7 = norm(vfinal)-sqrt((G*M)/rf);%sqrt((vxf^2)+(vyf^2))-sqrt((G*M)/rf);
b8 = dot(rfinal,vfinal)-0;
b9 = lambda_xf*vxf+lambda_yf*vyf - (G*M/(rf^3))*(lambda_vxf*xf+lambda_vyf*yf);

res = [b1;b2;b3;b4;b5;b6;b7;b8;b9];
return;

%------------------------------------------------------------%