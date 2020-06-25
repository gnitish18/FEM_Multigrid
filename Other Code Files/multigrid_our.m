function U = multigrid_our(Ah,f,v1,v2)
    [n1,~] = size(Ah);
    v(1:floor(n1/2),1) = 0;
    Ih(1:n1,1:floor(n1/2)) = 0;
    for i = 1:floor(n1/2)
        for j = 2*i-1:2*i+1
            if mod(j,2) == 0
                Ih(j,i) = 2;
            else
                Ih(j,i) = 1;
            end
        end
    end
    R = 0.5*Ih';
    A2h = project(Ah);
    f2h = R*f;
    V = Gauss_Siedel(A2h,f2h,v,v1);
%     [n2,~] = size(A2h)
%     v(1:floor(n2/2),1) = 0
%     I2h(1:n2,1:floor(n2/2)) = 0;
%     for i = 1:floor(n2/2)
%         for j = 2*i-1:2*i+1
%             if mod(j,2) == 0
%                 I2h(j,i) = 2;
%             else
%                 I2h(j,i) = 1;
%             end
%         end
%     end
%     R = 0.5*I2h';
%     A4h = project(A2h);
%     f4h = R*f2h ;
%     W = Gauss_Siedel(A4h,f4h,v,v1);
%     V = I2h*W;
      U = Ih*V;
end