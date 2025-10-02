clear
%forcing function
f = @(x) -exp(-x.^2/0.1);

%material properties
kappa = @(x) 0.75*sign(x)+0.85;

%boundary conditions
a = 0;
b = 0;

%domain size
xa = -1;
xb = 1;

L = xb - xa;  %length of the domain


%% We will solve d/dx (kappa(x) du/dx) = f(x) on domain x = [x_a,x_b] 
% with BC: 
%   u(x_a) = a
%   u(x_b)  = b
% using a central finite difference method. kappa(x) is a function of
% space.


N = 80; %number of internal gridpoints
h = L/(N+1); %interval length
M=N+2; %account for all points including boundary

%create grid 
x_grid = linspace(xa,xb,N+2); %this mesh includes endpoints
x = x_grid; 

%creade half-grid for kappa evaluations
x_kappa = (x(1:end-1) + x(2:end))/2;

%build matrix A
A = zeros(M,M); %initialize zero matrix NxN - this step is just for illustration, not necessary here

A(1,1) = h^2; %Dirichlet BC
for j=2:M-1
    A(j,j-1) = kappa(x_kappa(j-1));
    A(j,j) = -(kappa(x_kappa(j-1))+kappa(x_kappa(j)));
    A(j,j+1) = kappa(x_kappa(j));
end
A(M,M) = h^2; %Dirichlet BC
A = A/h^2

%NOTE: for large matrices it is a better idea to use sparse matrix
%representation

%create RHS vector
g = (f(x))';

%implement BC
g(1) = a;
g(M) = b;

%solve linear system
u = A\g;

%plot solution
xx= linspace(xa,xb);
plot(xx,-f(xx),'--k','LineWidth',2)
hold on;
plot(xx,kappa(xx),':k','LineWidth',2)
plot(x,u,'o','LineWidth',2,'MarkerSize',10);


hold off
xlabel('x'); ylabel('u(x)');
grid on
set(gca,'FontSize',18)

legend('forcing','kappa','solution')
