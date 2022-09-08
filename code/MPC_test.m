clear; close all; clc;
set(groot,'defaulttextinterpreter','Latex');
set(groot, 'defaultAxesTickLabelInterpreter','Latex');
set(groot, 'defaultLegendInterpreter','Latex');

%% Problem & Dynamics Setup

Var.lstar = 385692.5;           % characteristic length, km
Var.tstar = 377084.1526670386;  % Characteristic time, s

Var.m1 = 398600.4328969393;     % Gravitational parameter of Earth
Var.m2 = 4902.800582147765;     % Gravitational parameter of Luna
Var.rEL = 384400;   % Distance between Earth and Luna, nondim

Var.R1 = 6378;                  % Radius of Earth, km
Var.R2 = 1737;                  % Radius of Luna, km

Var.n = sqrt(Var.m1 / Var.rEL^3); % Mean motion of Luna
Var.D1 = Var.m2 * Var.rEL / (Var.m1 + Var.m2);
Var.mu = Var.m2 / (Var.m1 + Var.m2);
Var.L2x = L2LagrangePoint(Var.mu);

addpath('casadi-windows-matlabR2016a-v3.5.5');
import casadi.*
load('Shaloorbits.mat');

orbitNum = 10; % which orbit to use, 1 being closest to Lyapunov orbit

Traj_ref = x_Halo2{orbitNum};
Rx = Traj_ref(:,1); Ry = Traj_ref(:,2); Rz = Traj_ref(:,3);
Vx = Traj_ref(:,4); Vy = Traj_ref(:,5); Vz = Traj_ref(:,6);

numberOfOrbits = 5; % how many orbits to run optimal control

orbitalElementStep = 64; 
xd = [Rx'; Ry'; Rz'; Vx'; Vy'; Vz'];
w = waitbar(0, "preparing...");
t = Traj_ref(end,:);
%t = t(1:orbitalElementStep:end);
%xd = xd(1:orbitalElementStep:end, :)';

% Physical Constants
d_ES = 384400; % km, Earth-Sun distance
r_E = 6371; % km, radius of Earth
r_S = 695700; % km, radius of Sun
m1 = 4902.800582147765; % kg, mass of Earth
m2 = 398600.4328969393; % kg, mass of Sun
m_moon = 7.342e22; % kg, mass of Moon
mu = m1 / (m1 + m2); % mass ratio
%Omega = 1.99e-7; % rad/s, orbital angular velocity
%G = 6.67430e-20; % km^3*/kg/s^2;
%t_u = 1/Omega; % normalized time unit

sim_tim = t(end)*numberOfOrbits; % max simulation time
N = 32; % prediction horizon steps
window_time = t(end)/(size(xd,2)/N); % prediction horizon time
T = window_time/N; % sampling time
umax = 1; % maximum input magnitude

xd = repmat(xd, 1, numberOfOrbits+1); % repeats the desired state matrix for number of orbits

% Define the states
% SX.sym is casADi's version of symbolic variables
x = SX.sym('x'); 
y = SX.sym('y');
z = SX.sym('z');
vx = SX.sym('vx');
vy = SX.sym('vy');
vz = SX.sym('vz');
states = [x;y;z;vx;vy;vz];
n_states = length(states);

% Define the controls
ux = SX.sym('ux');
uy = SX.sym('uy');
uz = SX.sym('uz');
controls = [ux,uy,uz];
n_controls = length(controls);
P = SX.sym('P', n_states, N+3); % Parameters of initial and desired states 

% Equations of Motion (Transition Equations)
d = sqrt((x+mu)^2 + y^2 + z^2); % spacecraft distance from Sun
r = sqrt((x-1+mu)^2 + y^2 + z^2); % spacecraft distance from Earth
ax = 2*vy + x - (1-mu)*(x+mu)/(d^3) ...
    - mu*(x-1+mu)/(r^3) + ux; % acceleration in x
ay = -2*vx + y - (1-mu)*y/(d^3) ...
    - mu*y/(r^3) + uy; % acceleration in y
az = -(1-mu)*z/(d^3) - mu*z/(r^3) + uz; % acceleration in z
rhs = [vx; vy; vz; ax; ay; az]; % transition function

f = Function('f', {states,controls}, {rhs});
U = SX.sym('U', n_controls, N); % complete control matrix

X = SX.sym('X', n_states, (N+1)); % complete state matrix

X(:,1) = P(:,1); % set initial state

%% State Integration using fourth order Runge-Kutta
% calculates states from current state numerically

for k = 1:N
    st = X(:,k); % current state
    con = U(:,k); % current control
    if k == 1
        k1 = f(st, con + P(4:6, N+3));   % new 
        k2 = f(st + T/2*k1, con + P(4:6, N+3)); % new
        k3 = f(st + T/2*k2, con + P(4:6, N+3)); % new
        k4 = f(st + T*k3, con + P(4:6, N+3)); % new
    else
        k1 = f(st, con);   % new 
        k2 = f(st + T/2*k1, con); % new
        k3 = f(st + T/2*k2, con); % new
        k4 = f(st + T*k3, con); % new
    end
    st_next = st +T/6*(k1 + 2*k2 + 2*k3 + k4); % new state
    X(:,k+1) = st_next;
end

%% MPC Setup


% Complete transition function, yields all X given U
ff = Function('ff', {U,P}, {X});

% Initialize objective and constraints
obj = 0;
g = [];

% Define Q and R
Q = zeros(6,6); % Setup control bounds using Q. Leave all Q states
% as zero to provide an uncontrolled simulation.
Q(1,1) = 1; % x position
Q(2,2) = 1; % y position
Q(3,3) = 1; % z position
Q(4,4) = 0.01; % x velocity
Q(5,5) = 0.01; % y velocity
Q(6,6) = 0.01; % z velocity
Q = 1000*Q;

R = zeros(3,3);
R(1,1) = 1; % x control input
R(2,2) = 1; % y control input
R(3,3) = 1; % z control input
R = 1*R;

% Compute objective (cost) function
for k = 1:N
    st = X(:,k);
    con = U(:,k);
    obj = obj + (st-P(:,k+1))'*Q*(st-P(:,k+1)) + con'*R*con; % add stage cost
end
st = X(:,N+1);
obj = obj + (st-P(:,N+1))'*Q*(st-P(:,N+1)); % add terminal cost

% Compute constraints
for k = 1:N+1
    if k <= N
        g = [g; U(1,k)^2+U(2,k)^2+U(3,k)^2 - umax^2];
    else
        g = [g; 0];
    end
end

OPT_variables = reshape(U,n_controls*N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

% Solver options
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-7;
opts.ipopt.acceptable_obj_change_tol = 1e-7;
solver = nlpsol('solver', 'ipopt', nlp_prob, opts); % ipopt is interior point

args = struct; 
% inequality constraints
args.lbg = -inf;
args.ubg = inf;
% input constraints, these are set to zero because the inequality
% constraints already constrain the input
args.lbx = -inf;
args.ubx = inf;


%% MPC Initialization

t0 = 0; % initial time
x0 = xd(:,1);
%x0(1) = x0(1) + 0.00001;
xs = [0 0 0 0 0 0];
xx(:,1) = x0; % xx is the history of states
t(1) = t0; % set initial time vector to t0
u0 = zeros(N,3); % control input guess initialized to zero
mpciter = 0; % iteration number of the mpc
xx1 = []; % temporary storage variable
u_cl = []; % temporary storage variable

%% Simulation


while(mpciter < sim_tim/T)
    p = [zeros(6,1)];
    p2 = [zeros(3,1) ; 0.0001*randn(); 0.0001*randn(); 0.0001*randn()]; %
    
    % Uncomment for Gaussian noise perturbation
    %p2 = [zeros(6,1)];
    
    x0 = x0 + p2*T;
    args.p = [x0 xd(:,mpciter+1:mpciter+N+1) p]; % set the values of the parameter vector
    args.x0 = reshape(u0', 3*N, 1); % initial value of the optimization variables
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, 'lbg', ...
    args.lbg, 'ubg', args.ubg, 'p', args.p); 
    u = reshape(full(sol.x)', 3, N)';
    ff_value = ff(u', args.p);
    xx1(:,1:6,mpciter+1) = full(ff_value)';
    u_cl = [u_cl; u(1,:)];

    t(mpciter+1) = t0;
    [t0,x0,u0] = shift(T, t0, x0, 0*u, f);
    xx(:,mpciter+2) = x0;
    mpciter = mpciter + 1;
    waitbar(mpciter / (sim_tim/T), w, sprintf("calculating optimal control: %.2f%%", 100*mpciter / (sim_tim/T)));
    
    % If the spacecraft diverges off track too much, stops simulation
    if norm(x0(1:3)) > 1e7/Var.lstar % 10^7 km dimensionalized
        disp("Result diverged. Simulation terminated.")
        break;
    end
end
close(w);

%plotMPC = 1;
%draw_video_CR3BP_3D(t,xx,xx1,u_cl,xs,N,xd, plotMPC)

%% Plots


figure(2);
x = xx(1,:) * Var.lstar;
y = xx(2,:) * Var.lstar;
z = xx(3,:) * Var.lstar;
%plot(L1(1), L1(2), 'r*');
PlotSys(Var)
hold on; grid on
plot3(Rx * Var.lstar, Ry * Var.lstar, Rz * Var.lstar, 'm-')
plot3(x, y, z, 'b-');
th = linspace(0,2*pi,1000);
axis equal;
legend('Luna', 'Earth', '$L_2$', 'Reference Trajectory', 'Spacecraft Trajectory', 'Interpreter', 'Latex');
xlabel('$x$ (km)', 'Interpreter', 'Latex')
ylabel('$y$ (km)', 'Interpreter', 'Latex')
zlabel('$z$ (km)', 'Interpreter', 'Latex')
xlim([0.8,1.2]*Var.lstar)
ylim([-0.2,0.2]*Var.lstar)
zlim([-0.2,0.2]*Var.lstar)

udim = Var.lstar / Var.tstar^2;

t = t*Var.tstar / 60 / 24; % hours

figure(3);
subplot(311);
stairs(t,u_cl(:,1) * udim,'b','linewidth',1.5); 
ylabel('$u_x (m/s^2)$');
grid on;
subplot(312);
stairs(t,u_cl(:,2) * udim,'b','linewidth',1.5);
ylabel('$u_y (m/s^2)$');
grid on;
subplot(313);
stairs(t,u_cl(:,3) * udim,'b','linewidth',1.5);
xlabel('t (hrs)');
ylabel('$u_z (m/s^2)$');
grid on;

vdim = Var.lstar/Var.tstar;

figure(4);
subplot(311);
stairs(t,xx(4,1:end-1) * vdim,'b','linewidth',1.5); 
ylabel('$v_x (m/s)$');
grid on;
subplot(312);
stairs(t,xx(5,1:end-1) * vdim,'b','linewidth',1.5);
ylabel('$v_y (m/s)$');
grid on;
subplot(313);
stairs(t,xx(6,1:end-1) * vdim,'b','linewidth',1.5);
xlabel('t (hrs)');
ylabel('$v_z (m/s)$');
grid on;

u_mag = sqrt(u_cl(:,1).^2+u_cl(:,2).^2+u_cl(:,3).^2);
u_tot = sum(u_mag)*T;

figure(5);
plot(t, u_mag,'b','linewidth',1.5); 
xlabel('t (hrs)');
ylabel('$||u|| (m/s^2)$');
grid on;

%{
figure(6);
subplot(321);
plot(t,xx(1,1:end-1) - xd(1,1:size(xx,2)-1) - xd(1,size(xx,2)-1)','b'); 
ylabel('$e_{x}$', 'Interpreter', 'Latex', 'FontSize', 16);
grid on;
subplot(323);
plot(t,xx(2,1:end-1) - xd(2,1:size(xx,2)-1) ,'b'); 
ylabel('$e_{y}$', 'Interpreter', 'Latex', 'FontSize', 16);
grid on;
subplot(325);
plot(t,xx(3,1:end-1) - xd(3,1:size(xx,2)-1) ,'b'); 
ylabel('$e_{z}$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('time', 'Interpreter', 'Latex', 'FontSize', 12);
grid on;
subplot(322);
plot(t,xx(4,1:end-1) - xd(4,1:size(xx,2)-1) ,'b'); 
ylabel('$e_{vx}$', 'Interpreter', 'Latex', 'FontSize', 16);
grid on;
subplot(324);
plot(t,xx(5,1:end-1) - xd(5,1:size(xx,2)-1) ,'b'); 
ylabel('$e_{vy}$', 'Interpreter', 'Latex', 'FontSize', 16);
grid on;
subplot(326);
plot(t,xx(6,1:end-1) - xd(6,1:size(xx,2)-1) ,'b'); 
ylabel('$e_{vz}$', 'Interpreter', 'Latex', 'FontSize', 16);

grid on;

grid on;
xlabel('time', 'Interpreter', 'Latex', 'FontSize', 12);
%}

% figure(6);
% subplot(411);
% stairs(t,xx(1,1:128) - xd(1,1:128),'b','linewidth',1.5); 
% ylabel('ex');
% grid on;
% subplot(412);
% stairs(t,xx(2,1:128) - xd(2,1:128),'b','linewidth',1.5);
% ylabel('ey');
% grid on;
