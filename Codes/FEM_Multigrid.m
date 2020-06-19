clear; 
close all; 
clc;

m = input("Enter the number of elements: ");    % no. of elements
n = m+1;                                        % no. of nodes
h = 5/m;                                        % length of each element, length of beam = 5 units
k = 1/h*[1, -1; -1, 1];                         % individual stiffness matrix
Ah = zeros(n,n);                                % deflection matrix
U = ones(n,1);
term3 = zeros(n,1);

% Define overall stiffness matrix
for i = 1:m
    Ah(i,i)     = Ah(i,i)     + k(1,1);
    Ah(i,i+1)   = Ah(i,i+1)   + k(1,2);
    Ah(i+1,i)   = Ah(i+1,i)   + k(2,1);
    Ah(i+1,i+1) = Ah(i+1,i+1) + k(2,2);
end

% Term 3 part of F in AX=F

for i = 1:m
    j = i-1;
    t1 = -j^3/3 + j^2 - j + 1/h*((i^4-j^4)/12 - (i^3-j^3)/3 + (i^2-j^2)/2);
    t2 = i^3/3 - i^2 + i - 1/h*((i^4-j^4)/12 - (i^3-j^3)/3 + (i^2-j^2)/2);
    %t1 = ((3*j^4 + i^4)/12) - ((2*j^3 + i^3)/3) + ((j^2 + i^2)/2) - i*(j^3/3 - j^2 + j);
    %t2 = ((j^4 + 3*i^4)/12) - ((j^3 + 2*i^3)/3) + ((j^2 + i^2)/2) - j*(i^3/3 - i^2 + i);
    term3(i) = term3(i) - 150 * t1;
    term3(i+1) = term3(i+1) - 150 * t2;
end

F = zeros(n,1);                         % Force matrix
x = 5;                                  % Length of the beam
F(n) = 1*(-150*x + 150*x^2 - 50*x^3);  % Last term - slope
F = F + term3;
disp(F)
det(Ah)
v = input("Number of iterations: ");
U = multigrid_gs(Ah, F, v, v);
%U = Gauss_Siedel(Ah, F, U, v)
% U = Gauss_Siedel_inversion(Ah, F, U, v)

term3 = zeros(n,1);
term3(1) = -50 + 12.5/h*2^4 - 50/h*2^3 + 75*h;
term3(2) = 0 - 12.5/h*(2*2^4-3^4) + 50/h*(2*2^3-3^3) + 150*h;
term3(3) = 0 - 12.5/h*(2*3^4-4^4-2^4) + 50/h*(2*3^3-4^3-2^3) + 150*h;
term3(4) = 0 - 12.5/h*(2*4^4-5^4-3^4) + 50/h*(2*4^3-5^3-3^3) + 150*h;
term3(5) = 0 - 12.5/h*(2*5^4-6^4-4^4) + 50/h*(2*5^3-6^3-4^3) + 150*h;
term3(6) = -50*6^3 - 12.5/h*(2*6^4-7^4-5^4) + 50/h*(2*6^3-7^3-5^3) + 75*h + 150*6^2;
disp(term3);
% 
term3 = zeros(n,1);
for i = 1:m
    j = i-1;
    t2 = -50*i^3 + 12.5*j^4/h - i^4 + 150*i^2 - 50*j^3/h + 50*i^3/h +75*h;
    t1 = 50*j^3 - 12.5*j^4/h + i^4 - 150*j^2 + 50*j^3/h - 50*i^3/h +75*h;
    term3(i) = term3(i) + t1;
    term3(i+1) = term3(i+1) + t2;
end
disp(term3);
% disp(term3);

%Ah = [1,1,-2,1,3,-1;2,-1,1,2,1,-3;1,3,-3,-1,2,1;5,2,-1,-1,2,1;-3,-1,2,3,1,3;4,3,1,-6,-3,-2]
%F = [4;20;-15;-3;16;-27]
%Ah = [5,2,-1,-1,2,1;4,3,1,-6,-3,-2;1,1,-2,1,3,-1;2,-1,1,2,1,-3;1,3,-3,-1,2,1;-3,-1,2,3,1,3]
%F = [-3;-27;4;20;-15;16]
%Ah= [4,1,-1;2,7,1;1,-3,12]
%F = [3; 19; 31]