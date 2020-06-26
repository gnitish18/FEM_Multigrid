%% Function to solve the matrix equation A*x = B using iterative Gauss-Seidel relaxation method
%   Input Parameters:
%   * A  - Coefficient Matrix
%   * B  - Resultant Matrix
%   * X  - Variable Matrix
%   * fl - Flag variable for termination condition
%   * v  - Number of iterations
%   Output Parameters:
%   * X  - Solved Variable Matrix
%%

function X = Gauss_Seidel(A, B, X, fl, v)
    % Find the size of the matrix
    [n,~]    = size(A);
    
    % Iterate Gauss-Seidal for v iterations when flag is disabled
    if ~fl
        
        for k = 1:v
            % Store previous iteration values
            x_old = X;
            for i = 1:n
                sigma = 0;
                for j = 1:i-1
                    sigma = sigma + A(i,j)*X(j);
                end
                for j = i+1:n
                    sigma = sigma + A(i,j)*x_old(j);
                end
                X(i) = (1/A(i,i))*(B(i) - sigma);
            end
            
        end
    
    % Iterate Gauss-Seidal until convergence when flag is enabled
    else
        % Initial error
        normval = 1;
        
        % Iterate until error is greater than tolerance
        while normval>0
            % Store previous iteration values
            x_old = X;
            
            for i = 1:n
                sigma = 0;
                for j = 1:i-1
                    sigma = sigma + A(i,j)*X(j);
                end
                for j = i+1:n
                    sigma = sigma + A(i,j)*x_old(j);
                end
                X(i) = (1/A(i,i))*(B(i) - sigma);
            end
            % Update error and find its norm
            normval = norm(X-x_old);
        end
        
    end
    
end