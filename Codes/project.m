function b = project(a)
    [n,~] = size(a);
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
    R = 0.5*I';
    b = R*a*I;
end
    
