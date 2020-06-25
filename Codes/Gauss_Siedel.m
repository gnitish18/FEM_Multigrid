function x = Gauss_Siedel(A, b, x, flag, v)
    
    [n,~]    = size(A);                 % Find the size of the matrix
    
    % Iterate Gauss-Seidal for v iterations when flag is disabled
    if ~flag
        
        for k = 1:v
            x_old = x;                  % Store previous iteration values
            for i = 1:n
                sigma = 0;
                for j = 1:i-1
                    sigma = sigma + A(i,j)*x(j);
                end
                for j = i+1:n
                    sigma = sigma + A(i,j)*x_old(j);
                end
                x(i) = (1/A(i,i))*(b(i) - sigma);
            end
            
        end
    
    % Iterate Gauss-Seidal until convergence when flag is enabled
    else
        normval = 1;                    % Initial error
        
        % Iterate until error is greater than tolerance
        while normval>0
            x_old = x;                  % Store previous iteration values
            for i = 1:n
                sigma = 0;
                for j = 1:i-1
                    sigma = sigma + A(i,j)*x(j);
                end
                for j = i+1:n
                    sigma = sigma + A(i,j)*x_old(j);
                end
                x(i) = (1/A(i,i))*(b(i) - sigma);
            end
            normval = norm(x-x_old);    % Update error and find its norm
        end
        
    end
    
end