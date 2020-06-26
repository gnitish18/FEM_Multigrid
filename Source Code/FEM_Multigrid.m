%% Program to determine the deflection of a Cantilever-Beam under Uniform Distributed Load using Finite Element Method (FEM) and solved using Multigrid Gauss-Seidel method
%%

clear; 
close all; 
clc;

% Get the number of elements
disp("Grid independency is achieved at 30 elements");
disp("However, feel free to test with any number of elements");
m = input("Enter the number of elements: ");

% If the number of elements is odd, makes it even
% to satisfy multigrid transformation matrix sizes
m = m + mod(m,2);

fprintf("Initializing... ");
n     = m+1;                    % Initialize number of nodes
l     = 5;                      % Define length of the beam
h     = l/m;                    % Length of each element
k     = 1/h*[1, -1; -1, 1];     % Individual stiffness matrix
Ah    = zeros(n,n);             % Deflection matrix
F     = zeros(n,1);             % Force matrix

% Define overall stiffness matrix
for i = 1:m
    Ah(i,i)     = Ah(i,i)     + k(1,1);
    Ah(i,i+1)   = Ah(i,i+1)   + k(1,2);
    Ah(i+1,i)   = Ah(i+1,i)   + k(2,1);
    Ah(i+1,i+1) = Ah(i+1,i+1) + k(2,2);
end

% Force matrix (F) in A*X = F
for i = 1:m
    % Initialize the nodal locations
    xi = i*h; 
    xj = xi-h;
    
    % Initialize Common term
    t = (xi^4-xj^4)/12 - (xi^3-xj^3)/3 + (xi^2-xj^2)/2;
    % Initialize the terms
    T1 = -xj^3/3 + xj^2 - xj + 1/h*t;
    T2 = xi^3/3 - xi^2 + xi - 1/h*t;
    
    % Update the Force Matrix
    F(i)   = F(i)   + 150*T1;
    F(i+1) = F(i+1) + 150*T2;
end

% Add last term of force matrix - slope
F(n) = F(n) + (-150*l + 150*l^2 - 50*l^3);
fprintf(" Done");
pause(0.3)
fprintf(repmat('\b', 1, 21));

% Solve the equation using Multigrid Gauss-Seidel
U = Multigrid_TwoGrid(Ah, F);

% Display the results and analysis
Results(U);
