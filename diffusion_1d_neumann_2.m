clear
%forcing function
f = @(x) 0*x;

sigma = 0;
gamma = 1;
b = 2;
xa = 0;
xb = 1;

L = xb - xa;  %length of the domain

%exact solution
u_exact = @(x) exp(x) + (sigma-exp(xa))*(x - xb) + b - exp(xb);

%% We will solve d^2(u)/dx^2 = f(x) on domain x = [x_a,x_b] 
% with BC: 
%   u'(x_a) = sigma
%   u'(x_b)  = gamma
% using a central finite difference method

neumann_method= 3; %1 - first order one sided difference

N = 100; %number of internal gridpoints
h = L/(N+1); %interval length
M=N+2; %account for all points including boundary

%create grid 
x_grid = linspace(xa,xb,M); %this mesh includes endpoints
x = x_grid; 

%build matrix A
A = zeros(M,M); %initialize zero matrix N+2 x N+2 
A = diag(-2*ones(M,1))    + ...    % main diagonal
    diag(ones(M-1,1), -1) + ...    % upper diagonal
    diag(ones(M-1,1),  1);         % lower diagonal

%modify first row to represent BC

switch neumann_method
    case 1 
        A(1,1) = -h;
        A(1,2) = h;
        A(M,M-1) = h;
        A(M,M) = -h;
    case 2 
        A(1,1) = -h;
        A(1,2) = h;
        A(M,M-1) = h;
        A(M,M) = -h;
    case 3
        A(1,1) = 3*h/2;
        A(1,2) = -2*h;
        A(1,3) = h/2;
        A(M,M) = 3*h/2;
        A(M,M-1) = -2*h;
        A(M,M-2) = h/2;
end

A = A/h^2

%NOTE: for large matrices it is a better idea to use sparse matrix
%representation

%create RHS vector
g = (f(x))';

%implement BC
switch neumann_method
    case 1
        g(1) = sigma;
        g(M) = gamma;
    case 2
        g(1) = sigma + h/2*f(x(1));
        g(M) = gamma + h/2*f(x(M));
    case 3
        g(1) = sigma;
        g(M) = gamma;
end



g(M) = b;

%solve linear system
u = A\g;

%plot solution
xx= linspace(xa,xb);
%plot(xx,u_exact(xx),'--k','LineWidth',2)
%hold on;
plot(x,u,'o','LineWidth',2,'MarkerSize',10);
%plot(xa,a,'ob','LineWidth',2,'MarkerSize',10);
%plot(xb,b,'ob','LineWidth',2,'MarkerSize',10);

%hold off
xlabel('x'); ylabel('u(x)');
grid on
set(gca,'FontSize',18)

