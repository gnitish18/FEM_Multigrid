function U = Multigrid_TwoGrid(Ah, f)

    disp("Running Multigrid...");
    [n,~]    = size(Ah);
    u(1:n,1) = 0;
    U(1:n,1) = 0;
    v1       = 1;
    flag     = 0;
    I        = zeros(n,floor(n/2));

    while (U(1) < 0 || flag == 0)
        
        disp("Iterating Gauss-Seidel");
        flag = 1;
        v1   = v1 + 1;
        u1   = Gauss_Siedel(Ah, f, u, 0, v1);
        
        for i = 1:floor(n/2)
            for j = 2*i-1:2*i+1
                if mod(j,2) == 0
                    I(j,i) = 2;
                else
                    I(j,i) = 1;
                end
            end
        end
        
        R   = 0.5*I';
        rh  = f - Ah*u1;
        r2h = R*rh;
        A2h = Project_FineToCoarse(Ah);
        e2h = inv(A2h)*r2h;
        eh  = I*e2h;
        u1  = u1 - eh;
        U   = Gauss_Siedel(Ah, f, u1, 1, 1);
    end
    
end