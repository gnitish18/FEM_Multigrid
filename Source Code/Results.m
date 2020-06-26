%% Function to plot the Analytical solution and compare it with the FEM solution
%   Input Parameters: 
%   * U_fem - FEM solution of Deflection Matrix
%%

function Results(U_fem)

    fprintf("Displaying Results...")
    
    % Find the size of the fine stiffness matrix
    [n,~] = size(U_fem);
    
    % Generate points along the length of beam
    X = linspace(0,5,n);

    figure(1);      % Initialize the plot name
    plot(X, U_fem); % Plot the deflection of the beam vs its nodal location
    title('Plot of FEM Solution')
    ylabel('Deflection of Beam (units)') 
    xlabel('Length of Beam (units)')
    grid on

    % Analaytical solution of deflection
    U_ana = -75*X.^2 + 50*X.^3 - 12.5*X.^4;
    
    figure(2);      % Initialize the plot name
    plot(X, U_ana); % Plot the deflection of the beam vs its nadal location
    title('Plot of Analytical Solution')
    ylabel('Deflection of Beam (units)') 
    xlabel('Length of Beam (units)')
    grid on
    
    figure(3);      % Initialize the plot name
    plot(X, U_fem); % Plot the deflection of the beam vs its nadal location
    title('Analytical and FEM Comparison')
    ylabel('Deflection of Beam (units)') 
    xlabel('Length of Beam (units)') 
    hold on         % Multiple plots on same graph
    scatter(X, U_ana);
    legend('FEM','Analytical')
    hold off
    grid on
    
    fprintf(" Done");
    pause(0.3)
    fprintf(repmat('\b', 1, 26));
    
    % Display the animation of Beam Deflection
    Animate(U_fem, U_ana, X);
    
end 
