function x = Gauss_Siedel(A, b, x, flag, v)

    if flag == 0
        n = size(x,1);
        for k = 1:v
            x_old = x;
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
    else
        n = size(x,1);
        normval=1;
        while normval>0
            x_old = x;
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
            normval = norm(x-x_old);
        end
    end
end