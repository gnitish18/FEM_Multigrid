function U = Multigrid_TwoGrid(Ah, F)

    disp("Running Multigrid...");
    [n,~]    = size(Ah);                    % Find the size of the fine stiffness matrix
    u(1:n,1) = 0;                           % Initial DEflection matrix Initialization for Gauss-Seidal
    U(1:n,1) = 0;                           % Define Deflection matrix
    v1       = 1;                           % Number of iterations for Gauss-Seidal
    flag     = 1;                           % Flag to ensure first iteration
    I        = zeros(n,floor(n/2));         % Define Projection Operator matrix

    disp("Iterating Gauss-Seidel with optimum iterations for the given elements");
    % Loop to find the optimum number of iterations 
    % for application of initial Gauss-Seidel relaxation
    while (U(1) < 0 || flag == 1)
        
        flag = 0;                           % Disable flag
        v1   = v1 + 1;                      % Increment the number of iterations
        
        % Applying Gauss-Seidel relaxation method with v1 iterations
        U    = Gauss_Siedel(Ah, F, u, 0, v1);
        
        % Initialize the Projection Operator matrix as tri-diagonal
        for i = 1:floor(n/2)
            for j = 2*i-1:2*i+1
                if mod(j,2) == 0
                    I(j,i) = 2;
                else
                    I(j,i) = 1;
                end
            end
        end
        
        R   = 0.5*I';                       % Initialize Restriction Matrix
        rh  = F - Ah*U;                     % Compute Residue
        r2h = R*rh;                         % Project Residue from fine grid to coarse
        A2h = Project_FineToCoarse(Ah);     % Project Stiffness matrix from fine grid to coarse
        e2h = A2h\r2h;                      % Solve to find the error
        eh  = I*e2h;                        % Interpolate error from coarse grid to fine
        U  = U - eh;                        % Update the solution
        
        % Applying Gauss-Seidel relaxation method until convergence 
        % with updated solution of initial deflection matrix
        U   = Gauss_Siedel(Ah, F, U, 1, 1);
    end
    
end