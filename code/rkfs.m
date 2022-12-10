
function [y,errest,nsteps]=rkfs(f,t,y0,varargin)
% RKFS
%
% Implementation of various runge-kutta-methods for solving ordinary
% differencial equations with fixed step size.
%
% By default, the solver computes solutions only at the given points of
% time with no substeps. In this case, the integration steps exactly
% correspond to the time given time vector.
%
% Some methods yield an error estimate and are capable of adapting the step
% size, however this is done only by internally inserting substeps as
% required to meet the required tolerance. The solver tries to satisfy both
% AbsTol and RelTol. RelTol is compared to 2*|y1-y2] / (|y1| + |y2|), where
% y1 and y2 are the low-order and high-order solutions.
%
% For implicit methods, the nonlinear system is solved iteratively in a
% nested loop.
%
% Butcher tableaus are taken from:
% https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
%
% Input arguments:
%  f    - function (handle) to integrate, of the form [dy,Jdy] = f(t, y)
%         The jacobian Jdy is necessary for implicit methods. When
%         using explicit methods, the function is evaluated with one output
%         argument (dy) only.
%  t    - vector of time points to calculate the solutions
%  y0   - vector initial conditions
%
% Optional parameters:
%  'Method'   - specifies the runge-kutta method, see the switch-case list
%               default: 'dormand-prince'
%  'AbsTol'   - absolute tolerance, scalar or vector, for adaptive steps
%               default: inf
%  'RelTol'   - relative tolerance, scalar
%               default: 1e-2
%  'Adaptive' - whether or not to insert substeps to achieve the specified
%               accuracy
%               default: false
%  'MinStep'  - substeps are only inserted as long as the step is not
%               smaller than the specified value
%               default: 1e-4
%  'MaxIterations'
%             - maximum number of sub-iterations for solving the nonlinear
%               system (implicit methods only)
%               default: 10
%  'NonlinearTol'
%             - termination tolerance for nonlinear solver (implicit
%               methods only)
%               default: 1e-3
%
% Tillmann Stübler
% updated Jan 2021
% tillmann.stuebler@akka.eu
%
    function [yn,nsteps]=rkstep(t,y,dt)
        % perform runge-kutta step
        t
        if impl
            % for implicit methods, a nonlinear system must be solved.

            % first, make a reasonable initial guess for the slopes k_1,
            % k_2, ...
            K=cell(size(c));
            for j=1:numel(c)
                K{j}=f(t+c(j)*dt,y);
            end
            K=[K{:}];

            % now, iteratively solve the non-linear system
            for u=1:p.Results.MaxIterations
                k=cell(size(c));
                Jk=cell(size(c));

                % evaluate the derivatives and their jacobian
                for j=1:numel(c)
                    [k{j},Jk{j}]=f(t+c(j)*dt,y+(A(j,:)*K')'*dt);
                end

                % calculate a first-order correction term to improve the
                % solution
                M=dt*repmat(vertcat(Jk{:}),1,numel(c)).*repelem(A,n,n)-eye(numel(c)*n);
                deltaK=M\(K(:)-vertcat(k{:}));

                % assign the new approximate solution for slopes k_1, ...
                K(:)=K(:)+deltaK;

                if norm(deltaK)<p.Results.NonlinearTol*norm(K(:))
                    % terminate iteration
                    break
                end
            end

        else
            % perform function evaluations according to butcher tableau
            K=zeros(n,numel(c));
            for j=1:numel(c)
                K(:,j)=f(t+c(j)*dt,y+dt*K*A(j,:)');
            end
        end

        % compute the solution
        dy=K*b';
        yn=y+dt*dy(:,1);

        % the following is for error estimation and adaptive substeps only
        nsteps=1;
        

    end
p=inputParser;
p.FunctionName=mfilename;
p.addParameter('Method','dormand-prince');
p.addParameter('AbsTol',inf);
p.addParameter('RelTol',1e-2);
p.addParameter('Adaptive',false);
p.addParameter('MinStep',1e-4);
p.addParameter('MaxIterations',10); % sub iterations for implicit methods
p.addParameter('NonlinearTol',1e-3); % tolerance for sub iterations
p.parse(varargin{:});
% choose butler tableau for the selected method
switch p.Results.Method
    case 'euler' % euler forward
        A=0;
        b=1;

    case 'explicit-midpoint'
        A=[0 0;.5 0];
        b=[0 1];

    case 'heun'
        A=[0 0;1 0];
        b=[.5 .5];

    case 'ralston'
        A=[0 0;2/3 0];
        b=[.25 .75];

    case 'rk3'
        A=[0 0 0;.5 0 0;-1 2 0];
        b=[1/6 2/3 1/6];

    case 'heun3'
        A=[0 0 0;1/3 0 0;0 2/3 0];
        b=[.25 0 .75];

    case 'ssprk3'
        A=[0 0 0;1 0 0;.25 .25 0];
        b=[1/6 1/6 2/3];

    case 'rk4'
        A=[0 0 0 0;.5 0 0 0;0 .5 0 0;0 0 1 0];
        b=[1/6 1/3 1/3 1/6];

    case '38rule'
        A=[0 0 0 0;1/3 0 0 0;-1/3 0 0 0;1 -1 1 0];
        b=[1/8 3/8 3/8 1/8];

    case 'dormand-prince'
        A=[...
            0 0 0 0 0 0 0;...
            .2 0 0 0 0 0 0;...
            3/40 9/40 0 0 0 0 0;...
            44/45 -56/15 32/9 0 0 0 0;...
            19372/6561 -25360/2187 64448/6561 -212/729 0 0 0;...
            9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0;...
            35/384 0 500/1113 125/192 -2187/6784 11/84 0 ...
            ];
        b=[35/384 0 500/1113 125/192 -2187/6784 11/84 0;...
            5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];

    case 'rk-fehlberg'
        A=[...
            0 0 0 0 0 0;...
            .25 0 0 0 0 0;...
            3/32 9/32 0 0 0 0;...
            1932/2197 -7200/2197 7296/2197 0 0 0;...
            439/216 -8 3680/513 -845/4104 0 0;...
            -8/27 2 -3544/2565 1859/4104 -11/40 0 ...
            ];
        b=[16/135 0 6656/12825 28561/56430 -9/50 2/55;...
            25/216 0 1408/2565 2197/4104 -1/5 0];

    case 'heun-euler'
        A=[0 0;1 0];
        b=[.5 .5;1 0];

    case 'fehlberg'
        A=[0 0 0;.5 0 0;1/256 255/256 0];
        b=[1/256 255/256 0;...
            1/512 255/256 1/512];

    case 'bogacki-shampine'
        A=[0 0 0 0;.5 0 0 0;0 .75 0 0;2/9 1/3 4/9 0];
        b=[2/9 1/3 4/9 0;...
            7/24 1/4 1/3 1/8];

    case 'cash-karp'
        A=[0 0 0 0 0 0;...
            .2 0 0 0 0 0;...
            3/40 9/40 0 0 0 0;...
            3/10 -9/10 6/5 0 0 0;...
            -11/54 5/2 -70/27 35/27 0 0;...
            1631/55296 175/512 575/13824 44275/110592 253/4096 0];
        b=[37/378 0 250/621 125/594 0 512/1771;...
            2825/27648 0 18575/48384 13525/55296 277/14336 1/4];

    case 'backward-euler'
        A=1;
        b=1;

    case 'implicit-midpoint'
        A=.5;
        b=1;

    case 'crank-nicolson'
        A=[0 0;.5 .5];
        b=[.5 .5];

    case 'crouzeix'
        A=[.5+3^.5/6 0;-3^.5/3 .5+3^.5/6];
        b=[.5 .5];

    case 'qin-zhang'
        A=[.25 0;.5 .25];
        b=[.5 .5];

    case 'kraaijevanger-spijker'
        A=[.5 0;-.5 2];
        b=[-.5 1.5];

    case 'nørsett'
        % Nørsett's three-stage, 4th order diagonally implicit RK method
        x=1.0685790213016284067;
        A=[x 0 0;.5-x x 0;2*x 1-4*x x];
        b=[1 6*(1-2*x)^2-2 1]./(6*(1-2*x)^2);

    case 'dirk-3'
        % 4-stage, 3th-order diagonally implicit RK method
        A=[.5 0 0 0;1/6 .5 0 0;-.5 .5 .5 0;1.5 -1.5 .5 .5];
        b=[1.5 -1.5 .5 .5];

    case 'lobatto-iiia-2'
        A=[0 0;.5 .5];
        b=[.5 .5;1 0];

    case 'lobatto-iiia-4'
        A=[0 0 0;5/24 1/3 -1/24;1/6 2/3 1/6];
        b=[1/6 2/3 1/6;-.5 2 -.5];

    case 'lobatto-iiib-4'
        A=[1/6 -1/6 0;1/6 1/3 0;1/6 5/6 0];
        b=[1/6 2/3 1/6;-.5 2 -.5];

    case 'lobatto-iiic-2'
        A=[.5 -.5;.5 .5];
        b=[.5 .5;1 0];

    case 'lobatto-iiic-4'
        A=[1/6 -1/3 1/6;1/6 5/12 -1/12;1/6 2/3 1/6];
        b=[1/6 2/3 1/6;-.5 2 -.5];

    case 'lobatto-iiic*-2'
        A=[0 0;1 0];
        b=[.5 .5];

    case 'lobatto-iiic*-4'
        A=[0 0 0;.25 .25 0;0 1 0];
        b=[1/6 2/3 1/6];

    case 'radau-ia-3'
        A=[.25 -.25;.25 5/12];
        b=[.25 .75];

    case 'radau-ia-4'
        A=[1/9 (-1-6^.5)/18 (-1+6^.5)/18;1/9 11/45+7*6^.5/360 11/45-43*6^.5/360;1/9 11/45+43*6^.5/360 11/45-7*6^.5/360];
        b=[1/9 4/9+6^.5/36 4/9-6^.5/36];

    case 'radau-iia-3'
        A=[5/12 -1/12;3/4 .25];
        b=[.75 .25];

    case 'radau-iia-5'
        A=[11/45-7*6^.5/360 37/225-169*6^.5/1800 -2/225+6^.5/75;...
            37/225+169*6^.5/1800 11/45+7*6^.5/360 -2/225-6^.5/75;...
            4/9-6^.5/36 4/9+6^.5/36 1/9];
        b=[4/9-6^.5/36 4/9+6^.5/36 1/9];

    otherwise
        error('unsupported runge-kutta method: %s',p.Results.Method)

end
c=sum(A,2);
n=numel(y0);
dt=diff(t); % calculate step sizes
% initialize outputs
t=t(:);
y=nan(numel(t),numel(y0));
y(1,:)=y0;
nsteps=zeros(size(t));
errest=zeros(size(t));
impl=~all(~triu(A),'all'); % implicit solvers are strictly lower triangular
% finally, integrate the ode in a loop
for i=1:numel(t)-1
    [y(i+1,:),nsteps(i+1)]=rkstep(t(i),y(i,:)',dt(i));
end
end
