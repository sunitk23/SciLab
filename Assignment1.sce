disp("enter the matrix")
A = [input("value"), input("value"), input("value"); input("value"), input("value"), input("value"); input("value"), input("value"), input("value")];
disp("the matrix is:")
disp(A)
function Gauss(A)
    n= 3;
    for i=2:n
        for j=2:n
            A(i,j) = A(i,j) - A(1,j)*A(i,1)/A(1,1);
        end
        A(i,1) = 0;
    end
    
    for i=3:n
        for j=3:n
            A(i,j) = A(i,j) - A(2,j)*A(i,2)/A(2,2);
        end
        A(i,2) = 0;
    end
    disp("upper triangular matrix is")
    disp(A)
    disp(A(3,3), A(2,2), A(1,1), "the pivots are")
endfunction

function LU(A)
    U = A;
    m = det(U(1,1));
    n = det(U(2,1));
    a = n/m;
    U(2,:) = U(2,:) - U(1,:)/(m/n);
    n = det(U(3,1));
    b = n/m;
    U(3,:) = U(3,:) - U(1,:)/(m/n);
    m = det(U(2,2));
    n = det(U(3,2));
    c = n/m;
    U(3,:) = U(3,:) - U(2,:)/(m/n);
    disp(U, "the upper triangular matrix is:")
    L = [1,0,0;a,1,0;b,c,1];
    disp(L, "the lower triangular matrix is:")
endfunction

function Inverse(A)
    n = 3;
    Aug = [A, eye(n,n)]
    for  j=1:n-1
        for i = j+1:n
            Aug(i,j:2*n) = Aug(i,j:2*n) - Aug(i,j)/Aug(j,j)*Aug(j,j:2*n);
        end
    end
    for j=n: -1:2
        Aug(1:j-1,:) = Aug(1:j-1,:) - Aug(1:j-1,j)/Aug(j,j)*Aug(j,:);
    end
    for j = 1:n
        Aug(j,:) = Aug(j, :)/ Aug(j,j);
    end
    I = Aug(:,n+1:2*n);
    disp(I, "the inverse of the matrix is:")
endfunction





























