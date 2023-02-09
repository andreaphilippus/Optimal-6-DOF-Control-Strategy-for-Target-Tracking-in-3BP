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
Var.rEL = 384400;   % Distance between Earth and Moon, nondim

Var.R1 = 6378;                  % Radius of Earth, km
Var.R2 = 1737;                  % Radius of Moon, km

Var.n = sqrt(Var.m1 / Var.rEL^3); % Mean motion of Moon
Var.D1 = Var.m2 * Var.rEL / (Var.m1 + Var.m2);
Var.mu = Var.m2 / (Var.m1 + Var.m2);
Var.tol = 1e-12;                % Numerical Integration Tolerance

Var.NumOrb = 1;                 % Number of orbits to propagate

Var.vecSize = 5000;              % Attitude r   epresentation vector length, km

%% Obtain L2 Lagrange Point

Var.L2x = L2LagrangePoint(Var.mu);

%% Desired Condition
% Southern Halo
load('Shaloorbits.mat');

Traj_ref = x_Halo2{1500}(:,1:6);          % Reference Trajectory (pos/vel), nondim
tspan = x_Halo2{1500}(:,end);      % Reference time-span

%SB_ref = EA323toDCM(pi/6, -pi/3, pi/4);                 % Reference Attitude
SB_ref = eye(3);

%% Initial Condition

SB0 = SB_ref;
R0 = Traj_ref(1,1:3)';            % In N frame; with SB0 = eye(3), it is the same in B frame.
nu0 = Traj_ref(1,4:6)'; % In S frame
w0 = 0*1.3*[pi/4 pi/2 2*pi/6]';            % In B frame

for i = 1:length(Traj_ref)
    g_ref(:,:,i) = SE3(SB_ref, Traj_ref(i,1:3)');
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

%% MPC 

wb = waitbar(0, "preparing...");
tic

for i = 1:length(Traj_ref)
    xd(:,:,i) = SE3(SB_ref, Traj_ref(i,1:3)');
end

% Suppose that the Lunar Gateway does not rotate

vd = [zeros(3,length(Traj_ref)); Traj_ref(:,4:6)']; 

t = Traj_ref(end,:);
sim_tim = t(end) * Var.NumOrb;          % Maximum Simulation Time
N = 32;                                 % Prediction Horizon Step
window_time = t(end)/(size(xd,3)/N);    % Prediction Horizon Time
%T = window_time / N;   
T = linspace(0, window_time, N);                    % Sampling Time Steps (linear)

xd = cat(3, xd(:,:,1), repmat(xd(:,:,2:end), [1 1 Var.NumOrb+1]));      % Desired state, g^0
vd = [vd(:,1) repmat(vd(:,2:end), 1, Var.NumOrb+1)];                    % Desired velocity, v^0 

%% MPC Setup


% MPC Setup
t0 = 0;             % Initial time

%x0 = xd(:,:,1);     % Initial state
x0 = SE3(EA323toDCM(pi/4, pi/2, pi/5), xd(1:3,4,1) + 0*[0.002 0.001 -0.0005]');

%v0 = vd(:,1);       % Initial velocities
v0 = [w0; vd(4:6,1)];

xx(:,:,1) = x0;     % xx is the history of states
vv(:,1) = v0;       % vv is the history of velocities
t(1) = t0;          % Setting the initial time vector to t0
u0 = zeros(N,6);    % Six control inputs
mpciter = 0;        % Interation number of the mpc
xx1 = [];           % Temporary state variable storage
vv1 = [];           % Temporary velocity variable storage
u_cl = [];          % Temporary control variable storage
Q = 10000;
R = diag([1,1,1,10,10,10]);
timestep = T(2) - T(1);
% MPC
load('lindyn.mat')
load('Bf.mat')
while(mpciter < sim_tim/timestep)
    
    % perturbation

    rot_pert = 0.001 * randn(3,1);
    trs_pert = 0.005 * randn(3,1);

    v0 = v0 + [rot_pert; trs_pert] * timestep;

    % desired
    g_des(:,:,1) = x0;
    v_des(:,1) = v0;

    for i = mpciter+1:mpciter+N+6
        g_des(:,:,i-mpciter+1) = xd(:,:,i+1);
        v_des(:,i-mpciter+1) = vd(:,i+1);
    end
        
    % obtain A and B
    [c, r] = invSE3(x0);
    A = Afunc(c(1,1), c(1,2), c(1,3),...
        c(2,1), c(2,2), c(2,3),...
        c(3,1), c(3,2), c(3,3), ...
        Var.J(1,1), Var.J(2,2), Var.J(3,3),...
        Var.mu, v0(1), v0(2), v0(3),...
        r(1), r(2), r(3));
    B = Bfunc();

    cvx_begin quiet

        xtemp = [reshape(c,[],1); r; v0(4:6); v0(1:3)];
 
    % introduce variables:
        variable u(6,N-1);
        
        % cost function
        obj = 0;
        
        for i = 1:N-1
            BS_temp = reshape(xtemp(1:9,end), 3, 3);
            BS_err = g_des(1:3,1:3,i)' * BS_temp;

            R_temp = xtemp(10:12,end);
            R_err = R_temp - g_des(1:3,4,i);
            h_temp = SE3(BS_err, R_err);

            obj = obj + Q * power(2,norm(h_temp - eye(4), 'fro')) + u(:,i)'*R*u(:,i);
            
            xtemp = [xtemp, xtemp(:,i) + (A * xtemp(:,i) + B * u(:,i)) * timestep];
              
        end
        % final cost
        BS_temp = reshape(xtemp(1:9,end), 3, 3);
        BS_err = g_des(1:3,1:3,N)' * BS_temp;

        R_temp = xtemp(10:12,end);
        R_err = R_temp - g_des(1:3,4,N);
        h_temp = SE3(BS_err, R_err);
        obj = obj + Q * power(2,norm(h_temp - eye(4), 'fro'));
    
        minimize(obj)
        subject to 
            norm(u(:,1:3)) <= 0.5
            norm(u(:,4:6)) <= 10
        
    cvx_end
    
    %[ff_value1, ff_value2] = ff(u', args.p);
    
    %xx1(:,1:12,mpciter+1) = full(ff_value1)';   % Predicting Window
    %vv1(:,1:6,mpciter+1) = full(ff_value2)';    %          "      
    u_cl = [u_cl u(:,1)];
 
    %t(mpciter+1) = t0; % since linear, this is unnecessary
    
    %[t0,x0,v0,u0] = shiftCL(T, t0, LinST(x0), v0, u, f);
    
    [x0, v0] = shiftnonlinCVX(LinST(x0), v0, u(:,1), timestep, Var);

    x0 = SE3(reshape(x0(1:9),3,3), x0(10:12));

    xx(:,:,mpciter+2) = x0; % note that mpciter starts from 0
    vv(:,mpciter+2) = v0;

    
    mpciter = mpciter + 1;
    waitbar(mpciter / (sim_tim/timestep), wb, ...
        sprintf("calculating optimal control: %.2f%%", ...
        100*mpciter / (sim_tim/timestep)));
    
    
end
close(wb);
elaptime = toc;

%% MPC Result Analysis

for i = 1:size(xx,3)
    [SB_mpc(:,:,i), R_mpc(:,i)] = invSE3(xx(:,:,i));
end
R_mpc = R_mpc * Var.lstar;

t = linspace(0,sim_tim,size(R_mpc,2));
U_tau = u_cl(1:3,:); % Need dimensionalization quantity
U_f = u_cl(3:6,:) * Var.lstar / Var.tstar;

% Figureworks

%draw_video_bas(g_axisym0,g_ref0,Var,2)

%{
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
%title({['Trajectory of the Spacecraft under SE(3) Dynamics'],['Without Control Input'],['Target: Northern Halo Orbit']})
xlim([0.8,1.2]*Var.lstar)
ylim([-0.2,0.2]*Var.lstar)
zlim([-0.2,0.2]*Var.lstar)
%}
%draw_video_bas(g_axisym,g_ref,Var,4)

% MPC Figureworks
figure(8)

PlotSys(Var)

% Reference orbit
plot3(Traj_ref(:,1)*Var.lstar, Traj_ref(:,2)*Var.lstar, Traj_ref(:,3)*Var.lstar, 'b-')
plot3(R_mpc(1,:), R_mpc(2,:), R_mpc(3,:), 'r-')

for i = linspace(1, length(R_mpc), 8)
    DrawAttitude(floor(i), xx, Var, 'cl')
end

legend('Moon', 'Earth', '$L_2$',...
    'Reference Trajectory', 'Spacecraft Trajectory', 'location', 'best')

axis equal
xlabel('$\hat x$ (km)')
ylabel('$\hat y$ (km)')
zlabel('$\hat z$ (km)')
%title({['Trajectory of the Spacecraft under SE(3) Dynamics'],['Target: Northern Halo Orbit']})
xlim([0.8,1.2]*Var.lstar)
ylim([-0.2,0.2]*Var.lstar)
zlim([-0.2,0.2]*Var.lstar)

u_cl = u_cl';
tau_mag = sqrt(u_cl(:,1).^2 + u_cl(:,2).^2 + u_cl(:,3).^2);
f_mag = sqrt(u_cl(:,4).^2 + u_cl(:,5).^2 + u_cl(:,6).^2);

%%
figure(12)

for i = 1:size(xx,3)
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
t = t(1:end-1);

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

%% Video
%draw_video_CasADi(xx,xx1,N,xd,Var,7)
