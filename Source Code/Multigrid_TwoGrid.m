%% Function to solve the matrix equation using Multigrid Algorithm implementing two-grid
%   Input Parameters:
%   * Ah - Coefficient Fine Grid Matrix
%   * F  - Resultant Force Matrix
%   Output Parameters:
%   * U  - Solved Deflection Matrix
%%

function U = Multigrid_TwoGrid(Ah, F)

    fprintf("Running Multigrid...");
    % Find the size of the fine stiffness matrix
    [n,~]    = size(Ah); 
    
    u(1:n,1) = 0;       % Initialization of Initial Deflection matrix
    U(1:n,1) = 0;       % Define Final Deflection matrix
    U_0      = 0;       % Boundary condition
    v1       = 1;       % Number of iterations for Gauss-Seidal
    fl       = 1;       % Flag to ensure first iteration
    
    % Define Projection Operator matrix
    I = zeros(n,floor(n/2));         

    fprintf("\nIterating Gauss-Seidel with optimum iterations...");
    % Loop to find the optimum number of iterations 
    % for application of initial Gauss-Seidel relaxation
    % Error is U_0 - U(1), which is minimized
    while (U_0 - U(1) > 0 || fl == 1)
        
        fl = 0;         % Disable flag
        v1 = v1 + 1;    % Increment the number of iterations
        
        % Applying Gauss-Seidel relaxation method with v1 iterations
        U  = Gauss_Seidel(Ah, F, u, 0, v1);
        
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
        
        R   = 0.5*I';   % Initialize Restriction Matrix
        rh  = F - Ah*U; % Compute Residue
        r2h = R*rh;     % Project Residue from fine grid to coarse
        
        % Project Stiffness matrix from fine grid to coarse
        A2h = Project_FineToCoarse(Ah);
        
        e2h = A2h\r2h;  % Solve to find the error
        eh  = I*e2h;    % Interpolate error from coarse grid to fine
        U  = U - eh;    % Update the solution
        
        % Applying Gauss-Seidel relaxation method until convergence 
        % with updated solution of initial deflection matrix
        U   = Gauss_Seidel(Ah, F, U, 1, 1);
        
    end
    
    fprintf(" Done");
    pause(0.3)
    fprintf(repmat('\b', 1, 55));
    fprintf(" Done");
    pause(0.3)
    fprintf(repmat('\b', 1, 25));
    
end