clear; 
close all; 
clc;

% Get the number of elements
m = input("Enter the number of elements (even): ");

% If the number of elements is odd, makes it even
% to satisfy multigrid transformation matrix sizes
if mod(m,2)
    m=m+1;
    string = sprintf('The number of elements has been increased to %d as your value is not even', m);
    disp(string)
end

disp("Initializing... ");
n     = m+1;                    % Number of nodes
x     = 5;                      % Length of the beam
h     = x/m;                    % Length of each element
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
    xi = i*h;
    xj = xi-h;
    t1 = -xj^3/3 + xj^2 - xj + 1/h*((xi^4-xj^4)/12 - (xi^3-xj^3)/3 + (xi^2-xj^2)/2);
    t2 = xi^3/3 - xi^2 + xi - 1/h*((xi^4-xj^4)/12 - (xi^3-xj^3)/3 + (xi^2-xj^2)/2);
    F(i)   = F(i)   + 150*t1;
    F(i+1) = F(i+1) + 150*t2;
end

F(n) = F(n) + (-150*x + 150*x^2 - 50*x^3);          % Add last term of force matrix - slope
U = Multigrid_TwoGrid(Ah, F);                       % Solve the equation using Multigrid Gauss-Seidel

disp("Displaying Results...")
figure(1);                                          % Initialize the plot name
X = (0:m)*h;                                        % Define the x coordinates in terms of length of the beam
plot(X, U);                                         % Plot the deflection of the beam
