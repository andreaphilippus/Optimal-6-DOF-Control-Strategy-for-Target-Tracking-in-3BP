clear; close all; clc;
set(groot,'defaulttextinterpreter','Latex');
set(groot, 'defaultAxesTickLabelInterpreter','Latex');
set(groot, 'defaultLegendInterpreter','Latex');

%% Description
%{ 
This code produces Lyapunov and Halo orbit family about Earth-Luna L2.

Dependencies:
    CR3BP_EoM.m     <-- Equations of motion of CR3BP
    STMdot.m        <-- Gives time derivative of state transition matrix
    Diffcor.m       <-- Performs differential correction


Note that all values are nondimensionalized, except for m1 and m2.

Constant values are achieved from a dissertation by Pavlak (2013),
 "Trajectory Design and Orbit Maintenance Strategies in Multi-Body
 Dynamical Regimes".
%}

%% Variables

Var.m1 = 398600.4328969393;     % Gravitational parameter of Earth, kg
Var.m2 = 4902.800582147765;     % Gravitational parameter of Luna, kg

Var.R1 = 6378;                  % Radius of Earth, km
Var.R2 = 1737;                  % Radius of Luna, km
Var.lstar = 385692.5;           % characteristic length, km
Var.tstar = 377084.1526670386;  % Characteristic time, s
Var.mu = Var.m2 / (Var.m1 + Var.m2);
Var.tol = 1e-12;            % Numerical Integration Tolerance
Var.ds = 0.0005;             % Pseudo-arclength step size
Var.NumFam = 150;            % Number of orbits to generate for one family

%% Obtain L2 Lagrange Point

Var.L2x = L2LagrangePoint(Var.mu);

%% L1
%% Generating family of Northern Halo orbit
Var.ds = 0.0001;  
Var.NumFam = 600;

clc;

q0 = [0.8824; 0; -0.0618; 0; -0.2342; 0]; % Initial pos/vel guess
STM0 = eye(length(q0)); % Initial STM

for i = 1:Var.NumFam
    %fprintf('\nOrbit #: %d\n\n', i)
    Xguess = [q0; reshape(STM0, [], 1)]; % Initial state guess

    % Differential Correction
    [OutH.t(i), OutH.X(:,i), OutH.DF(:,:,i), OutH.MonoMat(:,:,i)] = Diffcor(Xguess, Var); 

    %disp(OutH.X(:,i))
    [t_Halo,s_Halo] = NumSolve(@(t,X)CR3BP_EoM(t,X,Var.mu), OutH.X(1:6,i), [0 2*OutH.t(i)], Var.tol, 0);
    
    x_Halo1N{i} = [s_Halo t_Halo];

    % Update x0, and use ydot0 value from the previous for better
    % differential correction itertaion process
    q0(1) = q0(1) - Var.ds;
    q0(3) = OutH.X(3,i);
    q0(5) = OutH.X(5,i);
end

%save('L1NorthHalo.mat', "x_Halo1N")

%% Generating family of Southern Halo orbit
clc;

q0 = [0.8824; 0; 0.0618; 0; -0.2342; 0]; % Initial pos/vel guess

for i = 1:Var.NumFam
    %fprintf('\nOrbit #: %d\n\n', i)
    Xguess = [q0; reshape(STM0, [], 1)]; % Initial state guess

    % Differential Correction
    [OutH.t(i), OutH.X(:,i), OutH.DF(:,:,i), OutH.MonoMat(:,:,i)] = Diffcor(Xguess, Var); 

    %disp(OutH.X(:,i))
    [t_Halo,s_Halo] = NumSolve(@(t,X)CR3BP_EoM(t,X,Var.mu), OutH.X(1:6,i), [0 2*OutH.t(i)], Var.tol, 0);
    
    x_Halo1S{i} = [s_Halo t_Halo];

    % Update x0, and use ydot0 value from the previous for better
    % differential correction itertaion process
    q0(1) = q0(1) - Var.ds;
    q0(3) = OutH.X(3,i);
    q0(5) = OutH.X(5,i);
end

%save('L1SouthHalo.mat', "x_Halo1S")


%% L2
%% Generating family of Lyapunov orbit

q0 = [Var.L2x + 0.01; 0; 0; 0; -0.005; 0]; % Initial pos/vel guess
STM0 = eye(length(q0)); % Initial STM

fo0"
    fprintf('\nOrbit #: %d\n\n', i)
    Xguess = [q0; reshape(STM0, [], 1)]; % Initial state guess

    % Differential Correction
    [Out.t(i), Out.X(:,i), Out.DF(:,:,i), Out.MonoMat(:,:,i)] = Diffcor(Xguess, Var); 
    [t_Lyap,s_Lyap] = NumSolve(@(t,X)CR3BP_EoM(t,X,Var.mu), Out.X(1:6,i), [0 2*Out.t(i)], Var.tol, 0);
    
    x_Lyap{i} = [s_Lyap t_Lyap];
    % Update x0, and use ydot0 value from the previous for better
    % differential correction itertaion process
    q0(1) = q0(1) + Var.ds;
    q0(5) = Out.X(5,i);
end

save('lyapunovorbits.mat', "x_Lyap")

%% Generating family of Northern Halo orbit
clc;

q0 = [1.18; 0; 0.0139; 0; -0.1; 0]; % Initial pos/vel guess

for i = 1:Var.NumFam
    fprintf('\nOrbit #: %d\n\n', i)
    Xguess = [q0; reshape(STM0, [], 1)]; % Initial state guess

    % Differential Correction
    [OutH.t(i), OutH.X(:,i), OutH.DF(:,:,i), OutH.MonoMat(:,:,i)] = Diffcor(Xguess, Var); 

    disp(OutH.X(:,i))
    [t_Halo,s_Halo] = NumSolve(@(t,X)CR3BP_EoM(t,X,Var.mu), OutH.X(1:6,i), [0 2*OutH.t(i)], Var.tol, 0);
    
    x_Halo{i} = [s_Halo t_Halo];

    % Update x0, and use ydot0 value from the previous for better
    % differential correction itertaion process
    q0(1) = q0(1) - Var.ds;
    q0(3) = OutH.X(3,i);
    q0(5) = OutH.X(5,i);
end

save('Nhaloorbits.mat', "x_Halo")

%% Generating family of Southern Halo orbit
clc;

q0 = [1.18; 0; -0.0139; 0; -0.1; 0]; % Initial pos/vel guess
%q0 = [1.140000000000000 0 -0.163379753195529 0 -0.223431558314127 0]';

for i = 1:Var.NumFam
    fprintf('\nOrbit #: %d\n\n', i)
    Xguess = [q0; reshape(STM0, [], 1)]; % Initial state guess

    % Differential Correction
    [OutH.t(i), OutH.X(:,i), OutH.DF(:,:,i), OutH.MonoMat(:,:,i)] = Diffcor(Xguess, Var); 

    disp(OutH.X(:,i))
    [t_Halo2,s_Halo2] = NumSolve(@(t,X)CR3BP_EoM(t,X,Var.mu), OutH.X(1:6,i), [0 2*OutH.t(i)], Var.tol, 0);
    
    x_Halo2{i} = [s_Halo2 t_Halo2];

    % Update x0, and use ydot0 value from the previous for better
    % differential correction itertaion process
    q0(1) = q0(1) - Var.ds;
    q0(3) = OutH.X(3,i);
    q0(5) = OutH.X(5,i);
end

save('Shaloorbits.mat', "x_Halo2")

%% Figurework
%% 

% Plot Lyapunov Orbit
figure(1)
PlotSys(Var)

for i = 1:10:length(x_Lyap)
    plot3(x_Lyap{i}(:,1), x_Lyap{i}(:,2), x_Lyap{i}(:,3), 'b-')
end

legend('Luna', 'Earth', '$L_2$', 'location', 'best')
axis equal
xlabel('$\hat x$')
ylabel('$\hat y$')
zlabel('$\hat z$')
xlim([0.4,1.6])
ylim([-0.6,0.6])

%Plot Northern Halo Orbit
figure(2)

PlotSys(Var)

for i = 1:50:length(x_Halo)
    plot3(x_Halo{i}(:,1), x_Halo{i}(:,2), x_Halo{i}(:,3), 'b-')
end
legend('Luna', 'Earth', '$L_2$', 'location', 'best')
axis equal
xlabel('$\hat x$')
ylabel('$\hat y$')
zlabel('$\hat z$')
xlim([0.75,1.25])
ylim([-0.25,0.25])
zlim([-0.25,0.25])


figure(3)

for i = 1:50:length(x_Halo1N)
    plot3(x_Halo1N{i}(:,1), x_Halo1N{i}(:,2), x_Halo1N{i}(:,3), 'b-')
end
axis equal
xlabel('$\hat x$')
ylabel('$\hat y$')
zlabel('$\hat z$')
xlim([0.75,1.25])
ylim([-0.25,0.25])
zlim([-0.25,0.25])

% Plot Southern Halo Orbit
figure(4)

PlotSys(Var)

for i = 1:50:length(x_Halo2)
    plot3(x_Halo2{i}(:,1), x_Halo2{i}(:,2), x_Halo2{i}(:,3), 'b-')
end
legend('Luna', 'Earth', '$L_2$', 'location', 'best')
axis equal
xlabel('$\hat x$')
ylabel('$\hat y$')
zlabel('$\hat z$')
xlim([0.75,1.25])
ylim([-0.25,0.25])
zlim([-0.25,0.25])

% Plot periodic orbits altogether
figure(5)

PlotSys(Var)

for i = 1:5:length(x_Lyap)
    plot3(x_Lyap{i}(:,1), x_Lyap{i}(:,2), x_Lyap{i}(:,3), 'b-')
end
for i = 1:50:length(x_Halo)
    plot3(x_Halo{i}(:,1), x_Halo{i}(:,2), x_Halo{i}(:,3), 'r-')
    plot3(x_Halo2{i}(:,1), x_Halo2{i}(:,2), x_Halo2{i}(:,3), 'g-')
end
legend('Luna', 'Earth', '$L_2$', 'location', 'best')
axis equal
xlabel('$\hat x$')
ylabel('$\hat y$')
zlabel('$\hat z$')
xlim([0.5,1.5])
ylim([-0.5,0.5])
ylim([-0.5,0.5])
