function [tCor, xCor, DFCor, STMCor] = Diffcor(Xguess, Var, Out)

maxIter = 50; % The code may perform up to 20 iterations
% If it does not converge to a value less than e-12 after 20
% iterations, send error and ask for a more appropriate initial guess.

X0 = Xguess;
mu = Var.mu;
tol = Var.tol;

% Take previous parameters for continuation
if nargin == 3
    Xl = Out.X(:,end);
    DFl = Out.DF(:,:,end);
end

fprintf('\nFirst Initial Condition Guess:\n')
disp(X0(1:6)')
% Time span is arbitrarily set
tspan = [0 20]; % nondimensionalized

Iter = 1; % Iteration counter
FX = [1;1;1]; % Temporary FX with its norm greater than the tolerance
while norm(FX) > tol
    if Iter > maxIter
        error('Max Iteration Reached. Try different initial guess.')
        break;
    end
   
    fprintf("\nDifferential Correction Iteration #: %d\n", Iter)
                  
    % Note that L2 orbit initial points are further away from the primary
    % bodies in +x direction. Therefore, an event is needed for a half
    % segment propagation, with direction +1.

    [t1, x1] = NumSolve(@(t,X)StateSTMdot(t, X, mu), X0, tspan, tol, 1);
    
    FX = [x1(end,2);x1(end,4);x1(end,6)];

    if nargin == 3 % Continuation
        FX = [FX; (1)];
    end
    

    fprintf('\nF(X):\n')
    disp(FX)

    fprintf('\n||F(X)||:\n')
    disp(norm(FX))

    % STM
    STM = reshape(x1(end,7:end),6,6);

    fprintf('\nSTM Φ:\n')
    disp(STM)

    % Jacobian Matrix

    % [xdot ydot zdot xddot yddot zddot] at final time for Jacobian Matrix
    Xdotf = CR3BP_EoM([], x1(end,1:6), mu);

    DF = [STM(2,3) STM(2,5) Xdotf(2);
        STM(4,3) STM(4,5) Xdotf(4);
        STM(6,3) STM(6,5) Xdotf(6)]; 

    %% Get update
    
    % Update X initial guess using Newton
    X0_update = -DF\FX; % This is the change in free variables

    fprintf('\nX(j+1):\n')
    disp(X0_update)

    % Switch only free variables, but set x0 fixed
    X0(3) = X0(3) + X0_update(1);
    X0(5) = X0(5) + X0_update(2);
    
    if norm(FX) > tol
        fprintf('\nNext Guess:\n')
        disp(X0(1:6)')

        Iter = Iter + 1;
    else
        fprintf('\nDifferential Correction Success!')
        fprintf('\nTotal number of iteration: %d\n', Iter)
    end
end

xCor = X0(1:6);
tCor = t1(end); % Multiply by 2 to get full orbital period.
DFCor = DF;
STMCor = STM;

end