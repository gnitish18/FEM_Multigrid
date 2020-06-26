%% Function to project a given Fine-Grid matrix to Coarse-Grid
%   Input Parameters: 
%   * A - Fine Grid Matrix
%   Output Paramerter:
%   * B - Coarse Grid Matrix
%%

function B = Project_FineToCoarse(A)
    % Find the size of the fine stiffness matrix
    [n,~] = size(A);
    % Define the Projection Operator matrix
    I(1:n,1:floor(n/2)) = 0;
    
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
    
    R = 0.5*(I');   % Initialize Restriction matrix
    B = R*A*I;      % Project the fine stiffness matrix to coarse

end 
