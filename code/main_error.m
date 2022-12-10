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

% Southern Halo
load('Shaloorbits.mat');

Traj_ref0 = x_Halo2{20};          % Reference Trajectory (pos/vel), nondim
tspan0 = x_Halo2{20}(:,end);      % Reference time-span

% Repeat for Var.NumOrb
Traj_ref = Traj_ref0;
tspan = tspan0;

for i = 1:Var.NumOrb-1
    Traj_ref = [Traj_ref(1:end-1,:); Traj_ref0];
    tspan = [tspan(1:end-1); tspan(end)+tspan0];
end

SB_ref = EA323toDCM(pi/6, -pi/3, pi/4);       % Reference Attitude
w_ref = [0 0 0]';       % Reference angular velocity

%% Halo Orbit
%% Initial Condition

SB0 = SB_ref;
R0 = Traj_ref(1,1:3)';            % In N frame; with SB0 = eye(3), it is the same in B frame.
nu0 = Traj_ref(1,4:6)'; % In N frame
w0 = [0 0 0]';

for i = 1:length(Traj_ref)
    g_ref(:,:,i) = SE3(SB0, Traj_ref(i,1:3)');
    v_ref(:,i) = [w0; Traj_ref(i,4:6)'];
end

X0 = [R0; w0; nu0; reshape(SB0,[],1)];      % Combined State for numerical integration

% Axisymmetric, No Control

Var.J = diag([1,1,3/2]);
Var.K = JtoK(Var.J);
Var.I = [Var.J zeros(3); zeros(3) eye(3)];

%% Classical Dynamics

[t_axisym, g_axisym, v_axisym] = propagator(X0, tspan, Var);

for i = 1:length(t_axisym)
    [SB_axisym(:,:,i), R_axisym(:,i)] = invSE3(g_axisym(:,:,i));
    EP_axisym(:,i) = DCMtoEP(SB_axisym(:,:,i));
end

R_axisym = R_axisym * Var.lstar;

%% Error Dynamics

g_d = g_ref;
g_l = g_d;

v_d = v_ref;
v_l = v_d;

% To be fixed later so that initial error can be introduced.
eta0 = zeros(6,1);
g0 = SE3(SB0, R0);
v0 = [w0; nu0];
input0 = [eta0; reshape(g0, [], 1); v0];

global ti
ti = 1;

%opt = odeset('Reltol',Var.tol,'Abstol',Var.tol);
%[t, out] = ode45(@(t,input)error_dyn(t,input,g_d,v_d,Var), tspan, input0);
% Use RK4?
%[out, ~,~] = rkfs(@(t,input)error_dyn(t,input,g_d,v_d,Var), tspan, input0);
[tout, out] = odeK4(@(t,input)error_dyn(t,input,g_d,v_d,Var),tspan,input0);

eta(1:6,:) = out(1:6,:);
for i = 1:length(out)
    temp = out(7:22,i);
    g(:,:,i) = reshape(temp,4,4);
    [~, temp] = invSE3(g(:,:,i));
    R(:,i) = temp;
end


%% Figureworks
%{
figure(1)

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
%}

figure(2)

PlotSys(Var)

% Reference orbit
plot3(Traj_ref(:,1)*Var.lstar, Traj_ref(:,2)*Var.lstar, Traj_ref(:,3)*Var.lstar, 'b-')
plot3(R(1,:)*Var.lstar, R(2,:)*Var.lstar, R(3,:)*Var.lstar, 'r-')

for i = linspace(1, length(R), 2*Var.NumOrb)
    DrawAttitude(floor(i), g, Var, 'cl')
end

legend('Moon', 'Earth', '$L_2$',...
    'Reference Trajectory', 'Spacecraft Trajectory', 'location', 'best')

axis equal
xlabel('$\hat x$ (km)')
ylabel('$\hat y$ (km)')
zlabel('$\hat z$ (km)')
title({['Trajectory of the Spacecraft under SE(3) Error Dynamics'],['Without Control Input'],['Target: Northern Halo Orbit']})
xlim([0.8,1.2]*Var.lstar)
ylim([-0.2,0.2]*Var.lstar)
zlim([-0.2,0.2]*Var.lstar)

