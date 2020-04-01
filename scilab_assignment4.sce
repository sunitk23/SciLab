disp("enter the matrix A")
A = [input("value"), input("value"), input("value"); input("value"), input("value"), input("value"); input("value"), input("value"), input("value")];
disp("the independent vectors stored on the columns of A:")
disp(A)

//Gram-Schmidt Orthogonalization in R3
function GSO(A)
    [m,n] = size(A);
    for k=1:n
        V(:,k) = A(:,k);
        for j=1:k-1
            R(j,k) = V(:,j)'*A(:,k);
            V(:,k) = V(:,k)-R(j,k)*V(:,j);
        end
        R(k,k) = norm(V(:,k));
        V(:,k) = V(:,k)/R(k,k);
    end
    disp('The Orthogonal basis is:')
    disp(V, 'Q=');
endfunction

//Find Eigen values and Eigen vectors of the given matrix
function Eigen(A)
    lam = poly(0, 'lam')
    charMat = A - lam*eye(3,3)
    disp(charMat, 'The charateristic matrix is:')
    charPoly = poly(A, 'lam')
    disp(charPoly, 'The charactersitic polynomial is:')
    lam = spec(A)
    disp(lam, 'The eigen values of A are:')
    function[x,lam] = eigenvectors(A)
        [n,m] = size(A);
        lam = spec(A)';
        x = [];
        for k=1:3
            B = A-lam(k)*eye(3,3); //characteristic matrix
            C = B(1:n-1,1:n-1); //coefficient matrix for the reduced system
            b = -B(1:n-1, n); // RHS vector for the reduced system
            y = C\b; // solution for the reduced system
            y = [y;1]; //complete eigen vector
            y = y/norm(y); //make unit eigen vector
            x = [x y];
        end
    endfunction
    [x, lam] = eigenvectors(A)
    disp(x, 'The eigen vectors of A are:');
endfunction

//Find largest Eigen value of A using Rayleigh Power method
function Largest_Eigen(A)
    u0 = [1 1 1]';
    disp(u0, 'The initial vector is')
    v = A*u0;
    a = max(u0);
    disp(a, 'First approximation to eigen value is:');
    while abs(max(v)-a) > 0.002
        disp(v, 'current eigen vector is:')
        a = max(v);
        disp(a, 'current eigen value is')
        u0 = v/max(v);
        v = A*u0;
    end
    format('v', 4);
    disp(max(v), 'The largest Eigen value is:')
    format('v', 5)
    disp(u0, 'The corresponding Eigen vector is:')
endfunction










