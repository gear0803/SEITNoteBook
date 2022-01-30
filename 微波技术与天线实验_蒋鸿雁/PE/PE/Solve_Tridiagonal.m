function U = Solve_Tridiagonal(main,sup,sub,rhs)

[~,n] = size(main);
new_sup = zeros(1,n-1);
new_rhs = zeros(1,n);

new_sup(1,1)=sup(1,1)/main(1,1);
for index = 2:n-1
    new_sup(1,index)=sup(1,index)/(main(1,index) - new_sup(1,index-1)*sub(1,index-1) );
end

new_rhs(1,1)=rhs(1,1)/main(1,1);
for index = 2:n
    new_rhs(1,index)=(rhs(1,index) - new_rhs(1,index-1)*sub(1,index-1))/...
        (main(1,index) - new_sup(1,index-1)*sub(1,index-1) );
end

U(1,n)=new_rhs(1,n);
for index = n-1 : -1 : 1
    U(1,index) = new_rhs(1,index) - new_sup(1,index)*U(1,index+1);
end
    



