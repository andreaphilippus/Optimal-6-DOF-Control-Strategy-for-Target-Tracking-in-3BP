%THIS IS GOOD CODE

clc
clear
close all

addpath('casadi-windows-matlabR2016a-v3.5.5');
import casadi.*

load('Shaloorbits.mat');

orbitNum = 10; % which orbit to use, 1 being closest to Lyapunov orbit

Traj_ref = x_Halo2{orbitNum};
Rx = Traj_ref(:,1); Ry = Traj_ref(:,2); Rz = Traj_ref(:,3);
Vx = Traj_ref(:,4); Vy = Traj_ref(:,5); Vz = Traj_ref(:,6);

numberOfOrbits = 2; %how many orbits to run optimal control
orbitalElementStep = 64; 
xd = [Rx'; Ry'; Rz'; Vx'; Vy'; Vz'];
w = waitbar(0, "preparing...");
t = Traj_ref(end,:);
%t = t(1:orbitalElementStep:end);
%xd = xd(1:orbitalElementStep:end, :);

%physical constants
d_ES = 384400; % km, Earth-Sun distance
r_E = 6371; % km, radius of Earth
r_S = 695700; % km, radius of Sun
m1 = 4902.800582147765; % kg, mass of Earth
m2 = 398600.4328969393; % kg, mass of Sun
m_moon = 7.342e22; % kg, mass of Moon
mu = m1 / (m1 + m2); % mass ratio
%Omega = 1.99e-7; %rad/s, orbital angular velocity
G = 6.67430e-20; %N*m^2/kg^2;
%t_u = 1/Omega; %normalized time unit

%Lagrange point coordinates
L1 = [1 - nthroot(mu/3, 3), 0, 0]';
L2 = [1 + nthroot(mu/3, 3), 0, 0]';
L4 = [(1/2)*(m2-m1)/(m1+m2), sqrt(3)/2, 0]';
L5 = [(1/2)*(m2-m1)/(m1+m2), -sqrt(3)/2, 0]';

sim_tim = t(end)*numberOfOrbits; %max simulation time
N = 32; %prediction horizon steps
window_time = t(end)/(size(xd,2)/N); %prediction horizon time
T = window_time/N; %sampling time
umax = 1; %
a_orbit = 40000; %orbit height for LEO ICs, km
th_orbit = pi/2; %orbit angle for LEO ICs, km

xd = repmat(xd, 1, numberOfOrbits+1); %repeats the desired state matrix for number of orbits

%define the states
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
vx = SX.sym('vx');
vy = SX.sym('vy');
vz = SX.sym('vz');
states = [x;y;z;vx;vy;vz];
n_states = length(states);

%define the controls
ux = SX.sym('ux');
uy = SX.sym('uy');
uz = SX.sym('uz');
controls = [ux,uy,uz];
n_controls = length(controls);
P = SX.sym('P', n_states, N+3); %parameters of initial and desired states 

%transition equations
d = sqrt((x+mu)^2 + y^2 + z^2);
r = sqrt((x-1+mu)^2 + y^2 + z^2);
ax = 2*vy + x - (1-mu)*(x+mu)/(d^3) ...
    - mu*(x-1+mu)/(r^3) + ux;
ay = -2*vx + y - (1-mu)*y/(d^3) ...
    - mu*y/(r^3) + uy;
az = -(1-mu)*z/(d^3) + mu*z/(r^3) + uz;
rhs = [vx; vy; vz; ax; ay; az]; %transition function

f = Function('f', {states,controls}, {rhs});
U = SX.sym('U', n_controls, N); %complete control matrix

X = SX.sym('X', n_states, (N+1)); %complete state matrix

X(:,1) = P(:,1); %set initial state

%define next states in terms of current states and controls
%RK4 discretization
for k = 1:N
    st = X(:,k); %current state
    con = U(:,k); %current control
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
    %st_next(1) = 1/2*T^2*
    %st_next(1:3) = st(1:3) + st(4:6)*T + (1/2)*f_value*T^2; %predict position
    %st_next(4:6) = st(4:6) + f_value*T;
    X(:,k+1) = st_next;
end

%complete transition function, yields all X given U
ff = Function('ff', {U,P}, {X});

%initialize objective and constraints
obj = 0;
g = [];

%define Q and R
Q = zeros(6,6);
Q(1,1) = 1; %x
Q(2,2) = 1; %y
Q(3,3) = 1; %z
Q(4,4) = 0.01; %vx
Q(5,5) = 0.01; %vy
Q(6,6) = 0.01; %vz
Q = 1000*Q;

R = zeros(3,3);
R(1,1) = 1; %ux
R(2,2) = 1; %uy
R(3,3) = 1; %uz
R = 1*R;

%compute objective (cost) function
for k = 1:N
    st = X(:,k);
    con = U(:,k);
    obj = obj + (st-P(:,k+1))'*Q*(st-P(:,k+1)) + con'*R*con; %add stage cost
end
st = X(:,N+1);
obj = obj + (st-P(:,N+1))'*Q*(st-P(:,N+1)); %add terminal cost

%compute constraints
for k = 1:N+1
    if k <= N
        g = [g; U(1,k)^2+U(2,k)^2+U(3,k)^2 - umax^2];
    else
        g = [g; 0];
    end
end

OPT_variables = reshape(U,n_controls*N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

%solver options
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-7;
opts.ipopt.acceptable_obj_change_tol = 1e-7;
solver = nlpsol('solver', 'ipopt', nlp_prob, opts); %ipopt is interior point

args = struct;
%inequality constraints
args.lbg = -inf;
args.ubg = 0;
%input constraints
args.lbx = -inf;
args.ubx = inf;

%{
x_orbit = (r_E+a_orbit)/d_ES*cos(th_orbit) + (1-mu);
y_orbit = (r_E+a_orbit)/d_ES*sin(th_orbit);
vx_orbit = -sqrt(G*m1/(r_E+a_orbit))*t_u/d_ES*sin(th_orbit);
vy_orbit = sqrt(G*m1/(r_E+a_orbit))*t_u/d_ES*cos(th_orbit);
%}

%MPC setup
t0 = 0; %initial time
%x0 = [x_orbit, y_orbit, 0, vx_orbit, vy_orbit, 0]'; %initial state, LEO
x0 = xd(:,1);

%x0(1) = 1; x0(2) = 0; x0(3) = 0;
%x0(5) = 0; x0(4) = 0; x0(6) = 0;
%x0(3) = 0.002
%x0(4) = -0.02;
%x0(1) = (1-x0(1))/2 + x0(1);
%x0 = [x_orbit-0.005, y_orbit, 0, 0, 0, 0]'; %initial state, LEO, no velocity

xs = [L1(1) L2(2) 0 0 0 0];
xx(:,1) = x0; %xx is the history of states
t(1) = t0; %set initial time vector to t0
u0 = zeros(N,3); %three control inputs
mpciter = 0; %iteration number of the mpc
xx1 = []; %temporary storage variable
u_cl = []; %temporary storage variable



while(mpciter < sim_tim/T)
%while(norm((x0-xs),2) > 1e-3 && mpciter < 2)
    %args.p = xd(:,mpciter+1:mpciter+N); %set the values of the parameter vector
    %p = [zeros(3,1) ; 0.05*randn(3,1)]; %random gaussian perturbation
%     if (mod(mpciter,10) == 1)
%         p2 = [zeros(3,1) ; 0.01*randn()+0.01; .01*randn(); 0.01*randn()]; %random gaussian perturbation
%     else
%         p2 = [zeros(6,1)];
%     end
    %p2 = [zeros(3,1) ; 0.01*sin(mpciter*T / (27.3*60*60*24/t_u)); 0.01*cos(mpciter*T / (27.3*60*60*24/t_u)); 0.01*cos(mpciter*T / (27.3*60*60*24/t_u))]; %random gaussian perturbation
    p2 = [zeros(3,1) ; 0.01*randn()+0.01; 0.01*randn(); 0.01*randn()];
    %p2 = [zeros(6,1)];
    p = [zeros(6,1)];
    x0 = x0 + p2*T;
    args.p = [x0 xd(:,mpciter+1:mpciter+N+1) p]; %set the values of the parameter vector
    args.x0 = reshape(u0', 3*N, 1); %initial value of the optimization variables
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, 'lbg', ...
    args.lbg, 'ubg', args.ubg, 'p', args.p);
    u = reshape(full(sol.x)', 3, N)';
    ff_value = ff(u', args.p);
    xx1(:,1:6,mpciter+1) = full(ff_value)';
    u_cl = [u_cl; u(1,:)];

    t(mpciter+1) = t0;
    [t0,x0,u0] = shift(T, t0, x0, u, f);
    xx(:,mpciter+2) = x0;
    mpciter = mpciter + 1;
    waitbar(mpciter / (sim_tim/T), w, sprintf("calculating optimal control: %.2f%%", 100*mpciter / (sim_tim/T)));
end
close(w);

%draw_video_CR3BP_3D(t,xx,xx1,u_cl,xs,N,xd,1)

figure(2);
x = xx(1,:);
y = xx(2,:);
plot(L2(1), L2(2), 'r*');
hold on;
plot3(xx(1,:), xx(2,:), xx(3,:), 'b-');
plot(1-mu, 0, 'go', 'Markersize', 2);
axis equal;
legend('Sun-Earth L1 Point', 'trajectory', 'Earth', 'Interpreter', 'Latex');
xlabel('x', 'Interpreter', 'Latex');
ylabel('y', 'Interpreter', 'Latex');
zlabel('z', 'Interpreter', 'Latex');

figure(3);
subplot(311);
stairs(t,u_cl(:,1),'b','linewidth',1.5); 
ylabel('ux');
grid on;
subplot(312);
stairs(t,u_cl(:,2),'b','linewidth',1.5);
ylabel('uy');
grid on;
subplot(313);
stairs(t,u_cl(:,3),'b','linewidth',1.5);
xlabel('time');
ylabel('uz');
grid on;

figure(4);
subplot(311);
stairs(t,xx(4,1:end-1),'b','linewidth',1.5); 
ylabel('vx');
grid on;
subplot(312);
stairs(t,xx(5,1:end-1),'b','linewidth',1.5);
ylabel('vy');
grid on;
subplot(313);
stairs(t,xx(6,1:end-1),'b','linewidth',1.5);
xlabel('time');
ylabel('vy');
grid on;

u_mag = sqrt(u_cl(:,1).^2+u_cl(:,2).^2+u_cl(:,3).^2);
u_tot = sum(u_mag)*T;

figure(5);
plot(t, u_mag,'b','linewidth',1.5); 
xlabel('time');
ylabel('u');

% figure(6);
% subplot(411);
% stairs(t,xx(1,1:128) - xd(1,1:128),'b','linewidth',1.5); 
% ylabel('ex');
% grid on;
% subplot(412);
% stairs(t,xx(2,1:128) - xd(2,1:128),'b','linewidth',1.5);
% ylabel('ey');
% grid on;
