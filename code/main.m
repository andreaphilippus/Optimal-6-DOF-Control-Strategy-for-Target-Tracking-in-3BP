clear; close all; clc;
set(groot,'defaulttextinterpreter','Latex');
set(groot, 'defaultAxesTickLabelInterpreter','Latex');
set(groot, 'defaultLegendInterpreter','Latex');

%% Description
%{ 


Dependencies:
    SE3.m               <-- Converts attitude and position state into one matrix on SE(3)
    se3alg.m            <-- se(3) Lie algebra
    propagator.m        <-- Performs numerical integration of the dynamics for a given timespan
    invSE3.m            <-- Gives separate attitude and position state from a matrix on SE(3)
    gravforceN.m        <-- Calculates gravitational force in N frame
    Dyn_Comb.m          <-- Function of dynamics for both SE(3) state matrix and velocity vector
    CrossProd.m         <-- Gives cross-product operator matrix on SO(3)
    Coadj.m             <-- Gives coadjoint matrix of a 6x1 vector
    Adj.m               <-- Gives adjoint matrix of a 6x1 vector
    L2LagrangePoint.m   <-- Gives L2 Lagrange Point
    PlotSys.m           <-- Plots Earth, Moon, and L2 Lagrange Point
    CasADiDynCont.m     <-- Function of dynamics for SE(3) state, velocity, and control
    draw_video_CasADi.m <-- Animates the motion of spacecraft 
    
    
%}

%% Variables

Var.lstar = 385692.5;           % characteristic length, km
Var.tstar = 377084.1526670386;  % Characteristic time, s

Var.m1 = 398600.4328969393;     % Gravitational parameter of Earth
Var.m2 = 4902.800582147765;     % Gravitational parameter of Moon
Var.rEL = 384400;               % Distance between Earth and Moon

Var.R1 = 6378;                  % Radius of Earth, km
Var.R2 = 1737;                  % Radius of Moon, km

Var.n = sqrt(Var.m1 / Var.rEL^3); % Mean motion of Moon
Var.D1 = Var.m2 * Var.rEL / (Var.m1 + Var.m2);
Var.mu = Var.m2 / (Var.m1 + Var.m2);
Var.tol = 1e-12;                % Numerical Integration Tolerance

Var.NumOrb = 2;                 % Number of orbits to propagate

Var.vecSize = 5000;              % Attitude r   epresentation vector length, km

%% Obtain L2 Lagrange Point

Var.L2x = L2LagrangePoint(Var.mu);

%% Desired Condition

% Lyapunov
load('lyapunovorbits.mat');

Traj_ref0 = x_Lyap{50};          % Reference Trajectory (pos/vel), nondim
tspan0 = x_Lyap{50}(:,end);      % Reference time-span

% Southern Halo
load('Shaloorbits.mat');

Traj_ref = x_Halo2{10};          % Reference Trajectory (pos/vel), nondim
tspan = x_Halo2{10}(:,end);      % Reference time-span

SB_ref = EA323toDCM(pi/6, -pi/3, pi/4);                 % Reference Attitude

%% Lyapunov Orbit
%% Initial Condition

SB0 = SB_ref;
R0 = Traj_ref0(1,1:3)';         % In S frame
nu0 = Traj_ref0(1,4:6)';        % In S frame
w0 = [pi/4 pi/2 2*pi/6]';            % In B frame

for i = 1:length(Traj_ref)
    g_ref0(:,:,i) = SE3(SB0, Traj_ref0(i,1:3)');
end

X0 = [R0; w0; nu0; reshape(SB0,[],1)];      % Combined State for numerical integration

% Axisymmetric, No Control

Var.J = diag([1,1,3/2]);
Var.K = JtoK(Var.J);
Var.I = [Var.J zeros(3); zeros(3) eye(3)]; 

[t_axisym0, g_axisym0, v_axisym0] = propagator(X0, tspan0, Var);

for i = 1:length(t_axisym0)
    [SB_axisym0(:,:,i), R_axisym0(:,i)] = invSE3(g_axisym0(:,:,i));
    EP_axisym0(:,i) = DCMtoEP(SB_axisym0(:,:,i));
end

R_axisym0 = R_axisym0 * Var.lstar;

%% Halo Orbit
%% Initial Condition

SB0 = SB_ref;
R0 = Traj_ref(1,1:3)';            % In N frame; with SB0 = eye(3), it is the same in B frame.
nu0 = Traj_ref(1,4:6)'; % In N frame

for i = 1:length(Traj_ref)
    g_ref(:,:,i) = SE3(SB0, Traj_ref(i,1:3)');
end

X0 = [R0; w0; nu0; reshape(SB0,[],1)];      % Combined State for numerical integration

% Axisymmetric, No Control

Var.J = diag([1,1,3/2]);
Var.K = JtoK(Var.J);
Var.I = [Var.J zeros(3); zeros(3) eye(3)];

[t_axisym, g_axisym, v_axisym] = propagator(X0, tspan, Var);

for i = 1:length(t_axisym)
    [SB_axisym(:,:,i), R_axisym(:,i)] = invSE3(g_axisym(:,:,i));
    EP_axisym(:,i) = DCMtoEP(SB_axisym(:,:,i));
end

R_axisym = R_axisym * Var.lstar;

%% MPC on CasADi

% Spacecraft Specification (axisymmetric)
Var.J = diag([1,1,3/2]);
Var.K = JtoK(Var.J);
Var.I = [Var.J zeros(3); zeros(3) eye(3)];

% Desired attitude
%SB_ref = EA323toDCM(pi/4, pi/2, pi/5);
SB_ref = eye(3);

addpath('casadi-windows-matlabR2016a-v3.5.5');
import casadi.*

wb = waitbar(0, "preparing...");

orbitalElementStep = 64;    
for i = 1:length(Traj_ref)
    xd(:,:,i) = SE3(SB_ref, Traj_ref(i,1:3)');
end

% Suppose that the Lunar Gateway does not rotate
vd = [zeros(3,length(Traj_ref)); Traj_ref(:,4:6)']; 

t = Traj_ref(end,:);
sim_tim = t(end) * Var.NumOrb;          % Maximum Simulation Time
N = 64;                                 % Prediction Horizon Step
window_time = t(end)/(size(xd,3)/N);    % Prediction Horizon Time
T = window_time / N;                    % Sampling Time

xd = cat(3, xd(:,:,1), repmat(xd(:,:,2:end), [1 1 Var.NumOrb+1]));      % Desired state, g^0
vd = [vd(:,1) repmat(vd(:,2:end), 1, Var.NumOrb+1)];                    % Desired velocity, v^0 

% Define states using CasADi symbolic
x = SX.sym('x'); y = SX.sym('y'); z = SX.sym('z');
C11 = SX.sym('C11'); C12 = SX.sym('C12'); C13 = SX.sym('C13');
C21 = SX.sym('C21'); C22 = SX.sym('C22'); C23 = SX.sym('C23');
C31 = SX.sym('C31'); C32 = SX.sym('C32'); C33 = SX.sym('C33');
%states = SE3([C11 C12 C13; C21 C22 C23; C31 C32 C33], [x;y;z]);
states = [C11 C12 C13 C21 C22 C23 C31 C32 C33 x y z]';
n_states = length(states);

wx = SX.sym('wx'); wy = SX.sym('wy'); wz = SX.sym('wz');
vx = SX.sym('vx'); vy = SX.sym('vy'); vz = SX.sym('vz');
vels = [wx wy wz vx vy vz]';
n_vels = length(vels);

% Define control
tx = SX.sym('tx'); ty = SX.sym('ty'); tz = SX.sym('tz');
fx = SX.sym('fx'); fy = SX.sym('fy'); fz = SX.sym('fz');
controls = [tx ty tz fx fy fz]';
n_controls = length(controls);


% Parameter vectors
% Parameter vectors have dimensions of number of states times initial + horizon step

P = SX.sym('P', n_states, N+1); % Parameters of initial and desired states
D = SX.sym('P', n_vels, N+1);   %              "                 vels

% Dynamics with control
[statedot, cstatedot] = CasADiDynCont(states, vels, controls, Var);

%f = Function('f', {states,vels,controls}, {statedot});
%h = Function('h', {states,vels,controls}, {cstatedot});
f = Function('f', {states,vels,controls}, {statedot, cstatedot});

% Complete Control Matrix
% U is defined so that U = [u(0) ... u(N-1)]
% Control input over the horizon
U = SX.sym('U', n_controls, N);    

% X and V: predictions will be stored in these variables
X = SX.sym('X', n_states, (N+1));   % Complete State Matrix
V = SX.sym('V', n_vels, (N+1));     % Complete vel Matrix

X(:,1) = P(:,1);                    % Set initial state
V(:,1) = D(:,1);                    % Set initial vel

%% State Integration using Fourth-Order Runge-Kutta

%{
for k = 1:N
    st = X(:,k);    % current state
    cste = V(:,k);   % current vel
    con = U(:,k);   % current control
    if k == 1     % For initial, need to consider the initial velocity
        [k1,l1] = f(st, cste, con + D(1:6, N+1));                    % new 
        [k2,l2] = f(st + T/2*k1, cste + T/2*l1, con + D(1:6, N+1));  % new
        [k3,l3] = f(st + T/2*k2, cste + T/2*l2, con + D(1:6, N+1));  % new
        [k4,l4] = f(st + T*k3, cste + T*l3, con + D(1:6, N+3));      % new
    else                  
        [k1,l1] = f(st, cste, con);                   % new 
        [k2,l2] = f(st + T/2*k1, cste + T/2*l1, con); % new
        [k3,l3] = f(st + T/2*k2, cste + T/2*l2, con); % new
        [k4,l4] = f(st + T*k3, cste + T*l3, con);     % new
    end
    st_next = st + T/6*(k1 + 2*k2 + 2*k3 + k4);     % new state
    cst_next = cste + T/6*(l1 + 2*l2 + 2*l3 + l4);   % new vel
    X(:,k+1) = st_next;
    V(:,k+1) = cst_next;
end
%}

for k = 1:N
    st = X(:,k);    % current state
    cste = V(:,k);   % current vel
    con = U(:,k);   % current control

    [k1,l1] = f(st, cste, con);                   % new 
    [k2,l2] = f(st + T/2*k1, cste + T/2*l1, con); % new
    [k3,l3] = f(st + T/2*k2, cste + T/2*l2, con); % new
    [k4,l4] = f(st + T*k3, cste + T*l3, con);     % new
    
    st_next = st + T/6*(k1 + 2*k2 + 2*k3 + k4);     % new state
    cst_next = cste + T/6*(l1 + 2*l2 + 2*l3 + l4);   % new vel
    X(:,k+1) = st_next;
    V(:,k+1) = cst_next;
end

%% MPC Setup

% Complete Transition Function, Yielding all X and V given U
ff = Function('ff', {U,[P;D]}, {X,V});

% Initialize objective and constraints
obj = 0;    % Objective (cost) Function
cst = [];   % Constraints     

% Define Q and R

%Q = zeros(4);
%Q(1:3,1:3) = 0.001;             % Attitude
%Q(1:3, 4) = 0.0001;             % Position

<<<<<<< Updated upstream
Q = diag([0.1 0.1 0.1]);

R = zeros(6);
R(1,1) = 0.01; % Torque in x
R(2,2) = 0.01; % Torque in y
R(3,3) = 0.01; % Torque in z
R(4,4) = 1; % Thrust in x
R(5,5) = 1; % Thrust in y
R(6,6) = 1; % Thrust in z
R = 1 * R;
Termweight = 10;

%Q = diag([1000 1000 1000]);
Q = 100;

S = diag([10 10 10 100 100 100]);

R = zeros(6);
R(1,1) = 0.1; % Torque in x
R(2,2) = 0.1; % Torque in y
R(3,3) = 0.1; % Torque in z
R(4,4) = 1; % Thrust in x
R(5,5) = 1; % Thrust in y
R(6,6) = 1; % Thrust in z

Termweight = 5;

% Compute the cost function
for k = 1:N
    % Add stage cost
    st_c = X(:,k);
    st = SE3(reshape(st_c(1:9), 3, 3), st_c(10:12)); % State, g

    % Reminder: P is the reference states
    p_c = P(:,k);
    p_test = SE3(reshape(p_c(1:9), 3, 3), p_c(10:12)); % Desired state, g0

    %p_test = SE3(reshape(p_c(1:9), 3, 3), p_c(10:12)); % Desired state, g0
    
    % test: reference state:
    %p_c = xd(:,:,k+1);

    % Velocity
    vel = V(:,k);
    d_c = D(:,k);

    con = U(:,k);
    
    % Error in attitude
    SB_e = reshape(p_c(1:9),3,3)' * reshape(st_c(1:9),3,3);
    
    % Error in position
    R_e = p_c(10:12) - st_c(10:12);

    % Error in velocity
    vel_e = d_c - vel;

    % Error in SE3
    h = SE3(SB_e, R_e); 
    
    %obj = obj + (st-P(:,k+1))'*Q*(st-P(:,k+1)) + con'*R*con; 

    %obj = obj + Q*norm(st - p_c, 'fro')^2 + con' * R * con;
    

    %obj = obj + norm(Q.*(h - eye(4)), 'fro')^2 + con' * R * con;

    %obj = obj + Q*norm((h - eye(4)), 'fro')^2 + vel_e' * S * vel_e + con' * R * con;
    obj = obj + Q*norm((h - eye(4)), 'fro')^2 + con' * R * con;

    
    %obj = obj + norm(R_e) + norm(SB_e - eye(3), 'fro') + con' * R * con;

    %obj = obj + Q * norm(p_c\st, 'fro')^2 + con' * R * con;

    %obj = obj + (st_c - p_c )
    
    %obj = obj + R_e' * Q * R_e + con' * R * con;
    
    obj = obj + 0;

end
% Add terminal cost
st_c = X(:,N+1);
st = SE3(reshape(st_c(1:9), 3, 3), st_c(10:12));

p_c = P(:,N+1);
p_test = SE3(reshape(p_c(1:9), 3, 3), p_c(10:12));

vel = V(:,N+1);
d_c = D(:,N+1);

SB_e = reshape(p_c(1:9),3,3)' * reshape(st_c(1:9),3,3);
R_e = p_c(10:12) - st_c(10:12);
vel_e = d_c - vel;

h = SE3(SB_e, R_e);


%obj = obj + (st-P(:,N+1))'*Q*(st-P(:,N+1)); 

%obj = obj + Termweight * norm(Q.*(h - eye(4)), 'fro')^2;

%obj = obj + R_e' * R_e;

%obj = 1000;

%obj = obj + (st-P(:,N+1))'*Q*(st-P(:,N+1)); 

%obj = obj + norm(Q*(h - eye(4)), 'fro')^2;

obj = obj + Q*norm((h - eye(4)), 'fro')^2;

%X_f = [R_e; reshape(SB_e - eye(3), 9, 1)];
%obj = obj + X_f' * X_f;


%obj = obj + Termweight * (norm(R_e) + norm(SB_e - eye(3), 'fro'));


% Compute the constraints

Var.f1theta = pi/4; % f1 constraint angle
Var.f2phi = pi/4;   % f2 constraint angle
Var.f3eta = pi/4;   % f3 constraint angle
Var.f4alpha = pi/4; % f4 constraint angle

Var.fmax = 30;       % Maximum orbital control input
Var.taumax = 5;     % Maximum attitude control input

for k = 1:N % Control constraint
    g = SE3(reshape(X(1:9,k),3,3), X(10:12,k));
    g0 = SE3(reshape(P(1:9,k),3,3), P(10:12,k));
    
    %cst = [cst; (g(1:3,1:3)'*(g0(1:3,4)-g(1:3,4)))'*[0 1 0]' + ...
    %    norm(g0(:,4)-g(:,4)) * cos(Var.f1theta)];           % f1
    %cst = [cst; (g0(1:3,1:3)'*(g0(1:3,4)-g(1:3,4)))'*[0 1 0]' + ...
    %    norm(g0(:,4)-g(:,4)) * cos(Var.f2phi)];             % f2
    %cst = [cst; -U(6,k) + norm(U(4:6,k)) * cos(Var.f3eta)]; % f3
    %cst = [cst; (g(1:3,1:3)'*g(1:3,4))'*[0 0 1]' + ...
    %    norm(g(:,4)) * cos(Var.f4alpha)];                   % f4
    cst = [cst; U(4,k)^2+U(5,k)^2+U(6,k)^2 - Var.fmax^2];   % f5
    cst = [cst; U(1,k)^2+U(2,k)^2+U(3,k)^2 - Var.taumax^2];   % f6
    %cst = [cst; 1.5*Var.R1/Var.lstar - norm(g(1:3,4) - [Var.mu 0 0]')]; % f7
    %cst = [cst; 1.3*Var.R2/Var.lstar - norm(g(1:3,4) - [1-Var.mu 0 0]')]; % f8
    %cst = [cst; 0];
end


OPT_variables = reshape(U, n_controls*N, 1); % Variable to optimize: input
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', cst, 'p', [P;D]);

% Solver options
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.acceptable_tol = Var.tol;
opts.ipopt.acceptable_obj_change_tol = Var.tol;
solver = nlpsol('solver', 'ipopt', nlp_prob, opts); % ipopt is interior point

args = struct;

% Inequality constraints, -inf < g <= 0
args.lbg = -inf;
args.ubg = 0;

args.lbx = -inf;
args.ubx = inf;

%% MPC Setup
t0 = 0;             % Initial time

%x0 = xd(:,:,1);     % Initial state
x0 = SE3(EA323toDCM(pi/4, pi/2, pi/5), xd(1:3,4,1) + [0.002 0.001 -0.0005]');
%x0 = SE3(eye(3), [1.3 0.2 0.02]');     

%v0 = vd(:,1);       % Initial velocities
v0 = [w0; nu0];

xx(:,:,1) = x0;     % xx is the history of states
vv(:,1) = v0;       % vv is the history of velocities
t(1) = t0;          % Setting the initial time vector to t0
u0 = zeros(N,6);    % Six control inputs
mpciter = 0;        % Interation number of the mpc
xx1 = [];           % Temporary state variable storage
vv1 = [];           % Temporary velocity variable storage
u_cl = [];          % Temporary control variable storage

% MPC

while(mpciter < sim_tim/T)

    % Perturbations
    % Position can't be perturbed (must be continuous)
    %pos_pert = 0.01 * randn(3,1);
   
    % Attitude also can't be
    %att_pert = EA323toDCM(0.0001 * randn(3,1));

    % Velocities are the dynamical quantities that can be perturbed
    rot_pert = 0.0005 * randn(3,1) * timestep;
    trs_pert = 0.001 * randn(3,1) * timestep;

    %x0 = SE_addition(x0, SE3(att_pert, [0 0 0]'));
    v0 = v0 + [rot_pert; trs_pert];

    % Set the values of the parameter vector
    
    % Known P has a structure of [init, ref(1), ref(2), ..., ref(N)]
    
    % Station-keeping scenario
    
    args.p = LinST(x0);
    for i = mpciter+1:mpciter+N+1
        args.p(:,i-mpciter) = LinST(xd(:,:,i-mpciter+1));
    end
    
    % Orbit insertion scenario
    %{
    args.p = LinST(x0);
    for i = mpciter+1:mpciter+N+6
        args.p(:,i-mpciter+1) = LinST(xd(:,:,1));
    end
    %}

    % Likewise, for V

<<<<<<< Updated upstream
    args.p = [args.p; vd(:,mpciter+1:mpciter+N+1)];
=======
    % Station-keeping 
    args.p = [args.p; v0 vd(:,mpciter+2:mpciter+N+7)];
>>>>>>> Stashed changes

    % Orbit insertion
    %{
    pv = v0;
    for i = mpciter+1:mpciter+N+6
        pv(:,i-mpciter+1) = vd(:,1);
    end
    
    args.p = [args.p; pv];
    %}
    
    %args.p = [LinST(x0) LinST(xd(:,:,mpciter+1:mpciter+N+1)) p; ...
    %    v0 vd(:,mpciter+1:mpciter+N+1) p2];  

    args.x0 = reshape(u0, 6*N, 1);      % Initial value of the optimization variables
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, 'lbg', ...
        args.lbg, 'ubg', args.ubg, 'p', args.p);
    
    % Uncontrol
    %u = zeros(6,N)';
    % Control
    u = reshape(full(sol.x)', 6, N)';

    [ff_value1, ff_value2] = ff(u', args.p);

    xx1(:,1:12,mpciter+1) = full(ff_value1)';   % Predicting Window
    vv1(:,1:6,mpciter+1) = full(ff_value2)';    %          "      
    u_cl = [u_cl; u(1,:)];

    t(mpciter+1) = t0;
    [t0,x0,v0,u0] = shiftCL(T, t0, LinST(x0), v0, u, f);
    x0 = SE3(reshape(x0(1:9),3,3), x0(10:12));
    xx(:,:,mpciter+2) = x0;
    vv(:,mpciter+2) = v0;
    
    mpciter = mpciter + 1;
    waitbar(mpciter / (sim_tim/T), wb, ...
        sprintf("calculating optimal control: %.2f%%", ...
<<<<<<< Updated upstream
        100*mpciter / (sim_tim/T)));
=======
        100*mpciter / (sim_tim/T)), 'interpreter', 'latex');
    
    % If the spacecraft diverges off track too much, stops simulation
    if norm(x0(1:3,4) - xd(1:3,4,i-mpciter+1)) > 1e10/Var.lstar % 10^10 km dimensionalized
        disp("Result diverged. Simulation terminated.")
        break;
    end
    
>>>>>>> Stashed changes
end
close(wb);

%% MPC Result Analysis

for i = 1:size(xx,3)
    [SB_mpc(:,:,i), R_mpc(:,i)] = invSE3(xx(:,:,i));
    EP_mpc(:,i) = DCMtoEP(SB_mpc(:,:,i));
end

R_mpc = R_mpc * Var.lstar;

t = t * Var.tstar;
U_tau = u_cl(:,1:3); % Need dimensionalization quantity
U_f = u_cl(:,3:6) * Var.lstar / Var.tstar;

%% Figureworks

%{
figure(1)

PlotSys(Var)

% Reference orbit
plot3(Traj_ref0(:,1)*Var.lstar, Traj_ref0(:,2)*Var.lstar, Traj_ref0(:,3)*Var.lstar, 'b-')
plot3(R_axisym0(1,:), R_axisym0(2,:), R_axisym0(3,:), 'r-')

for i = linspace(1, length(R_axisym0), 2*Var.NumOrb)
    DrawAttitude(floor(i), g_axisym0, Var, 'cl')
end

legend('Moon', 'Earth', '$L_2$',...
    'Reference Trajectory', 'Spacecraft Trajectory', 'location', 'best')
axis equal
xlabel('$\hat x$ (km)')
ylabel('$\hat y$ (km)')
zlabel('$\hat z$ (km)')
%title({['Trajectory of the Spacecraft under SE(3) Dynamics'],['Without Control Input'],['Target: Lyapunov Orbit']})
xlim([0.8,1.2]*Var.lstar)
ylim([-0.2,0.2]*Var.lstar)
zlim([-0.2,0.2]*Var.lstar)

%draw_video_bas(g_axisym0,g_ref0,Var,2)

figure(3)

PlotSys(Var)

% Reference orbit
plot3(Traj_ref(:,1)*Var.lstar, Traj_ref(:,2)*Var.lstar, Traj_ref(:,3)*Var.lstar, 'b-')
plot3(R_axisym(1,:), R_axisym(2,:), R_axisym(3,:), 'r-')

for i = linspace(1, length(R_axisym), 2*Var.NumOrb)
    DrawAttitude(floor(i), g_axisym, Var, 'cl')
end

legend('Moon', 'Earth', '$L_2$',...
    'Reference Trajectory', 'Spacecraft Trajectory', 'location', 'best')

axis equal
xlabel('$\hat x$ (km)')
ylabel('$\hat y$ (km)')
zlabel('$\hat z$ (km)')
title({['Trajectory of the Spacecraft under SE(3) Dynamics'],['Without Control Input'],['Target: Northern Halo Orbit']})
xlim([0.8,1.2]*Var.lstar)
ylim([-0.2,0.2]*Var.lstar)
zlim([-0.2,0.2]*Var.lstar)

%draw_video_bas(g_axisym,g_ref,Var,4)

%{
figure(5)

PlotSys(Var)

% Reference orbit
plot3(Traj_ref(:,1)*Var.lstar, Traj_ref(:,2)*Var.lstar, Traj_ref(:,3)*Var.lstar, 'b-')
plot3(R_sph(1,:), R_sph(2,:), R_sph(3,:), 'r-')

for i = linspace(1, length(R_sph), 2*Var.NumOrb)
    DrawAttitude(floor(i), g_sph, Var, 'cl')
end

legend('Moon', 'Earth', '$L_2$',...
    'Reference Trajectory', 'Spacecraft Trajectory', 'location', 'best')

axis equal
xlabel('$\hat x$ (km)')
ylabel('$\hat y$ (km)')
zlabel('$\hat z$ (km)')
%title({['Trajectory of the Spacecraft under SE(3) Dynamics'],['Without Control Input'],['Target: Northern Halo Orbit']})
xlim([0.8,1.2]*Var.lstar)
ylim([-0.2,0.2]*Var.lstar)
zlim([-0.2,0.2]*Var.lstar)
%}
%}

% MPC Figureworks

%draw_video_CasADi(xx,xx1,N,xd,Var,7)

figure(8)

PlotSys(Var)

% Reference orbit
plot3(Traj_ref(:,1)*Var.lstar, Traj_ref(:,2)*Var.lstar, Traj_ref(:,3)*Var.lstar, 'k-')

% Uncontrolled
load('Uncontrol_MPC.mat')
plot3(R_uncontrolled(1,:), R_uncontrolled(2,:), R_uncontrolled(3,:), 'b--')

% MPC Controlled
plot3(R_mpc(1,:), R_mpc(2,:), R_mpc(3,:), 'r-')

% Attitudes
for i = linspace(1, length(R_uncontrolled), 5)
    DrawAttitude(floor(i), xx_uncontrolled, Var, 'bl')
end

for i = linspace(1, length(R_mpc), 15)
    DrawAttitude(floor(i), xx, Var, 'cl')
end

legend('Moon', 'Earth', '$L_2$',...
    'Reference Trajectory', 'Uncontrolled Behavior', 'Controlled Behavior', 'location', 'best')

axis equal
xlabel('$\hat x$ (km)')
ylabel('$\hat y$ (km)')
zlabel('$\hat z$ (km)')
%title({['Trajectory of the Spacecraft under SE(3) Dynamics'],['Target: Northern Halo Orbit']})
xlim([0.8,1.2]*Var.lstar)
ylim([-0.2,0.2]*Var.lstar)
zlim([-0.2,0.2]*Var.lstar)

figure(9)   % f Input

subplot(311);
stairs(t,u_cl(:,4),'b','linewidth',1.5); 
ylabel('$f_x$');
grid on;
subplot(312);
stairs(t,u_cl(:,5),'b','linewidth',1.5);
ylabel('$f_y$');
grid on;
subplot(313);
stairs(t,u_cl(:,6),'b','linewidth',1.5);
xlabel('time (s)');
ylabel('$f_z$');
grid on;

figure(10)  % tau Input

subplot(311);
stairs(t,u_cl(:,1),'b','linewidth',1.5); 
ylabel('$\tau_x$');
grid on;
subplot(312);
stairs(t,u_cl(:,2),'b','linewidth',1.5);
ylabel('$\tau_y$');
grid on;
subplot(313);
stairs(t,u_cl(:,3),'b','linewidth',1.5);
xlabel('time (s)');
ylabel('$\tau_z$');
grid on;

figure(11)

subplot(211)
tau_mag = sqrt(u_cl(:,1).^2 + u_cl(:,2).^2 + u_cl(:,3).^2);

plot(t, tau_mag, 'b', 'linewidth', 1.5); grid on
ylabel('$||\tau||$')

subplot(212)
f_mag = sqrt(u_cl(:,4).^2 + u_cl(:,5).^2 + u_cl(:,6).^2);

plot(t, f_mag, 'b', 'linewidth', 1.5); grid on
ylabel('$||f||$')

%sgtitle('Control Input Magnitude vs Time')

figure(12)

for i = 1:size(xx,3)-1
    temp(:,:,i) = h_calc(xx(:,:,i), xd(:,:,i));
    pos_err(:,i) = temp(1:3,4,i);
    abs_pos_err(i) = norm(pos_err(:,i));
    
    att_err(:,:,i) = temp(1:3,1:3,i);
    EA323_e(1:3,i) = DCMtoEA323(att_err(:,:,i));
    
    er(i) = norm(temp(:,:,i) - eye(4), 'fro')^2;
end

plot(t, er); grid
xlabel('time')
ylabel('$||h-I_4||_F$')
%title('Error ($||h-I_4||_F$) vs Time')

figure(13)
subplot(3,6,[1:3; 7:9; 13:15])
plot(t, abs_pos_err); grid
xlabel('time')
ylabel('$||R_e||$')
subplot(3,6,4:6)
plot(t, reshape(temp(1,4,:),length(t),1)); grid
ylabel('$x_e$')
subplot(3,6,10:12)
plot(t, reshape(temp(2,4,:),length(t),1)); grid
ylabel('$y_e$')
subplot(3,6,16:18)
plot(t, reshape(temp(3,4,:),length(t),1)); grid
ylabel('$z_e$')
xlabel('time')
%sgtitle('Errors vs Time')

figure(14)
subplot(6,6,[1:3 7:9 13:15])
plot(t, f_mag, 'b', 'linewidth', 1.5); grid         % Change to stairs?
ylabel('$||f||$')
subplot(6,6,4:6)
stairs(t,u_cl(:,4),'b','linewidth',1.5); grid
ylabel('$f_x$')
subplot(6,6,10:12)
stairs(t,u_cl(:,5),'b','linewidth',1.5); grid
ylabel('$f_y$')
subplot(6,6,16:18)
stairs(t,u_cl(:,6),'b','linewidth',1.5); grid
ylabel('$f_z$');
subplot(6,6,19:21)
stairs(t,u_cl(:,1),'b','linewidth',1.5); grid
ylabel('$\tau_x$')
subplot(6,6,25:27)
stairs(t,u_cl(:,2),'b','linewidth',1.5); grid 
ylabel('$\tau_y$')
subplot(6,6,31:33)
stairs(t,u_cl(:,3),'b','linewidth',1.5); grid
ylabel('$\tau_z$')
xlabel('time')
subplot(6,6,[22:24 28:30 34:36])
plot(t, tau_mag, 'b', 'linewidth', 1.5); grid         % Change to stairs?
ylabel('$||\tau||$')
xlabel('time')
%sgtitle('Control Input vs Time')
