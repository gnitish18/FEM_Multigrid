function b = Project_FineToCoarse(a)
    
    [n,~] = size(a);            % Find the size of the fine stiffness matrix
    I(1:n,1:floor(n/2)) = 0;    % Define the Projection Operator matrix
    
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
    
    R = 0.5*(I');               % Initialize Restriction matrix
    b = R*a*I;                  % Project the fine stiffness matrix to coarse

end
    
