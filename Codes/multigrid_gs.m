function U = multigrid_gs(Ah,f,v1,v2)
    [n,~] = size(Ah);
    u(1:n,1) = 0;
    u1 = Gauss_siedel(Ah,f,u,v1);
    rh = f-Ah*u1;
    r2h = project(rh);
    A2h = project(Ah);
    e2h = inv(A2h)*r2h;
    disp(e2h);
    e2h = r2h/A2h;
    disp(e2h);
    I(1:n,1:floor(n/2)) = 0;
    for i = 1:floor(n/2)
        for j = 2*i-1:2*i+1
            if mod(j,2) == 0
                I(j,i) = 2;
            else
                I(j,i) = 1;
            end
        end
    end
    eh = I*e2h;
    u1 = u1+eh;
    U = Gauss_siedel(Ah,f,u1,v2);
end