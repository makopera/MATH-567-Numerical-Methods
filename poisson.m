%%Poisson equation solver

%exact solution
u_exact = @(x,y) sin(pi*x).*sin(pi*y);
f = @(x,y) -2*pi^2*sin(pi*x).*sin(pi*y);



%plot the exact solution
% [X,Y] = meshgrid(0:0.01:1,0:0.01:1);
% Z = u_exact(X,Y);
% surf(X,Y,Z)

%2D to 1D numbering
intma = @(i,j,N) (N+2)*j + i;
N_array = [10,20,40,80,160]

for k=1:length(N_array)
    %create the matrix
    N = N_array(k);
    M = (N+2)^2;

    x = linspace(0,1,N+2);
    y = linspace(0,1,N+2);
    h = x(2)-x(1);

    A = sparse(M,M);
    g = zeros(M,1);
    spy(A)

    %loop over the interior points
    for i=1:N
        for j=1:N
            I = intma(i,j,N)+1;
            I_lft = intma(i-1,j,N)+1;
            I_rgt = intma(i+1,j,N)+1;
            I_bot = intma(i,j-1,N)+1;
            I_top = intma(i,j+1,N)+1;
            A(I,I_lft) = 1;
            A(I,I_rgt) = 1;
            A(I,I_bot) = 1;
            A(I,I_top) = 1;
            A(I,I) = -4;

            g(I) = f(x(i+1),y(j+1));
        end
    end

    %boundary conditions
    %bottom
    j=0;
    for i=0:N+1
        I = intma(i,j,N)+1;
        A(I,I) = h^2;
        g(I) = 0;
    end

    %top
    j=N+1;
    for i=0:N+1
        I = intma(i,j,N)+1;
        A(I,I) = h^2;
        g(I) = 0;
    end

    %left
    i=0;
    for j=0:N+1
        I = intma(i,j,N)+1;
        A(I,I) = h^2;
        g(I) = 0;
    end


    %right
    i=N+1;
    for j=0:N+1
        I = intma(i,j,N)+1;
        A(I,I) = h^2;
        g(I) = 0;
    end

    A = A/h^2;

    u = A\g;
    U = zeros(N+2,N+2);
    for i=0:N+1
        for j=0:N+1
            I = intma(i,j,N);
            U(i+1,j+1) = u(I+1);
        end
    end
    spy(A)

    uu_exact = zeros(M,1);
    for i=0:N+1
        for j=0:N+1
            I = intma(i,j,N);
            uu_exact(I+1) = u_exact(x(i+1),y(j+1));
        end
    end

    error(k) = norm(u-uu_exact)/norm(uu_exact)

end
%[XX,YY] = meshgrid(x,y);
%surf(XX,YY,U)

loglog(N_array,error,'o-')
grid on