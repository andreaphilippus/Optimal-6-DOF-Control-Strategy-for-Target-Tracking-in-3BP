clear; close all; clc;
set(groot,'defaulttextinterpreter','Latex');
set(groot, 'defaultAxesTickLabelInterpreter','Latex');
set(groot, 'defaultLegendInterpreter','Latex');

%% Description
%{ 


Dependencies:
    SE3.m           <-- Converts attitude and position state into one matrix on SE(3)
    se3alg.m        <-- se(3) Lie algebra
    propagator.m    <-- Performs numerical integration of the dynamics for a given timespan
    invSE3.m        <-- Gives separate attitude and position state from a matrix on SE(3)
    gravforceN.m    <-- Calculates gravitational force in N frame
    Dyn_Comb.m      <-- Function of dynamics for both SE(3) state matrix and velocity vector
    CrossProd.m     <-- Gives cross-product operator matrix on SO(3)
    Coadj.m         <-- Gives coadjoint matrix of a 6x1 vector
    Adj.m           <-- Gives adjoint matrix of a 6x1 vector
    
%}

%% Variables

Var.lstar = 385692.5;           % characteristic length, km
Var.tstar = 377084.1526670386;  % Characteristic time, s

Var.m1 = 398600.4328969393;     % Gravitational parameter of Earth
Var.m2 = 4902.800582147765;     % Gravitational parameter of Luna
Var.rEL = 384400;   % Distance between Earth and Luna, nondim

Var.R1 = 6378;                  % Radius of Earth, km
Var.R2 = 1737;                  % Radius of Luna, km

Var.n = sqrt(Var.m1 / Var.rEL^3); % Mean motion of Luna, rad/s

Var.D1 = Var.m2 * Var.rEL / (Var.m1 + Var.m2);
Var.mu = Var.m2 / (Var.m1 + Var.m2);
Var.tol = 1e-12;                % Numerical Integration Tolerance

Var.NumOrb = 1;                 % Number of orbits to propagate

%% Obtain L2 Lagrange Point

Var.L2x = L2LagrangePoint(Var.mu);

%% Desired Condition

% Southern Halo
load('Shaloorbits.mat');

Traj_ref = x_Halo2{10};          % Reference Trajectory (pos/vel), nondim
R_ref = Traj_ref(:,1:3)*Var.lstar; % Reference position, km
nu_ref = Traj_ref(:,4:6)*Var.lstar/Var.tstar; % Reference velocity, km/s
tspan = x_Halo2{10}(:,end)*Var.tstar;      % Reference time-span, s

for i = 1:Var.NumOrb
    tspan = [tspan; tspan(2:end) + tspan(end)];
    R_ref = [R_ref; R_ref(2:end,:)];
    nu_ref = [nu_ref; nu_ref(2:end,:)];
end

% With gamma0 = 0, R0 and nu0 are same in both inertial and rotating frame.
R0 = R_ref(1,:)';         
nu0 = nu_ref(1,:)';
%nu0(2) = -nu0(2);

w0 = 0*[-0.1 0.05 0.1]';                  % In B frame

% with OB0 = eye(3), R0 and nu0 are same in all inertial/rotating and body
% fixed frame

NB0 = eye(3);

% Spacecraft Body Shape Configuration (in Body fixed frame)
Var.J = diag([1,1,1]);
Var.I = [Var.J zeros(3); zeros(3) eye(3)];

%% Orbits in Inertial Frame

[R_E, R_L, R_L2] = SysInert(tspan, Var);

for i = 1:length(tspan)

    NO(:,:,i) = NODCM(tspan(i), Var);

    R_ref_N(i,:) = NO(:,:,i) * R_ref(i,:)'; % Reference position in Inertial frame, km
end


%% Spacecraft Dynamics in Inertial Frame

for i = 1:length(tspan)
    F_G_N(i,:) = GravF_Inert(R_ref_N(i,:), R_E(i,:), R_L(i,:), Var);
    F_G_O(i,:) = NO(:,:,i)' * F_G_N(i,:)';
end

%% Integrate

% Initial state
X0 = [R0; w0; nu0; reshape(NB0,[],1)];

opt = odeset('RelTol',Var.tol, 'AbsTol', Var.tol);

tint = linspace(0, Var.NumOrb* tspan(end), 2000);
 
[t,y] = ode45(@(t,X) Dynamics(t,X,Var), tint, X0, opt);

R = y(:,1:3);
for i = 1:length(tint)
    NO(:,:,i) = NODCM(tint(i), Var);
    R_O(i,:) = NO(:,:,i) * R(i,:)';
end



%% Figurework

figure(1)

plot(R_E(:,1), R_E(:,2), "Color",  [91, 207, 244] / 255); hold on; grid
plot(R_L(:,1), R_L(:,2), "Color",[0.7 0.7 0.7])
plot(R_L2(:,1), R_L2(:,2), 'r-')
plot3(R_ref_N(:,1), R_ref_N(:,2), R_ref_N(:,3))
plot3(R(:,1),R(:,2),R(:,3), 'b-')
axis equal
xlabel('$\hat X$')
ylabel('$\hat Y$')
zlabel('$\hat Z$')
legend('Earth', 'Luna', '$L_2$',...
    'Reference Trajectory', 'Spacecraft Trajectory', 'location', 'best')

figure(2)
plot(tspan, F_G_N); hold on;
plot(tspan, F_G_O)
legend('X','Y','Z','x','y','z')

figure(3)

PlotSys(Var)
plot3(R_ref(:,1), R_ref(:,2), R_ref(:,3))
plot3(R_O(:,1), R_O(:,2), R_O(:,3))

axis equal
xlabel('$\hat x$')
ylabel('$\hat y$')
zlabel('$\hat z$')
xlim([-2,2]*Var.lstar)
ylim([-2,2]*Var.lstar)
zlim([-2,2]*Var.lstar)