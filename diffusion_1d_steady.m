clear
%forcing function
f = @(x) exp(x);

a = 1;
b = 2;
xa = 0;
xb = 1;

L = xb - xa;  %length of the domain

%exact solution
C = ((b-a)-(exp(xb)-exp(xa)))/L; %integration constant
D = a - exp(xa) - C*xa;
u_exact = @(x) exp(x) + C*x + D;

%% We will solve d^2(u)/dx^2 = f(x) on domain x = [x_a,x_b] 
% with BC: 
%   u(x_a) = a
%   u(x_b)  = b
% using a central finite difference method


N = 5; %number of internal gridpoints
h = L/(N+1); %interval length

%create grid 
x_grid = linspace(xa,xb,N+2); %this mesh includes endpoints
x = x_grid(2:N+1); %exclude the endpoints here

%build matrix A
A = zeros(N,N); %initialize zero matrix NxN - this step is just for illustration, not necessary here
A = diag(-2*ones(N,1))    + ...    % main diagonal
    diag(ones(N-1,1), -1) + ...    % upper diagonal
    diag(ones(N-1,1),  1);         % lower diagonal

A = A/h^2

%NOTE: for large matrices it is a better idea to use sparse matrix
%representation

%create RHS vector
g = (f(x))';

%implement BC
g(1) = g(1) - a/h^2;
g(N) = g(N) - b/h^2;

%solve linear system
u = A\g;

%plot solution
xx= linspace(xa,xb);
plot(xx,u_exact(xx),'--k','LineWidth',2)
hold on;
plot(x,u,'o','LineWidth',2,'MarkerSize',10);
plot(xa,a,'ob','LineWidth',2,'MarkerSize',10);
plot(xb,b,'ob','LineWidth',2,'MarkerSize',10);

hold off
xlabel('x'); ylabel('u(x)');
grid on
set(gca,'FontSize',18)

%% Convergence study

N = 5;
n_tests = 5;
error = zeros(1,n_tests);
h_list = zeros(1,n_tests); 
for k=1:n_tests
    
    h = L/(N+1); %interval length
    
    x_grid = linspace(xa,xb,N+2); %this mesh includes endpoints
    x = x_grid(2:N+1); %exclude the endpoints here

    A = zeros(N,N); %initialize zero matrix NxN - this step is just for illustration, not necessary here
    A = diag(-2*ones(N,1))    + ...    % main diagonal
        diag(ones(N-1,1), -1) + ...    % upper diagonal
        diag(ones(N-1,1),  1);         % lower diagonal
    A = A/h^2;

    %create RHS vector
    g = (f(x))';

    %implement BC
    g(1) = g(1) - a/h^2;
    g(N) = g(N) - b/h^2;

    %solve linear system
    u = A\g;

    %calculate error
    h_list(k) = h;
    u_ref = u_exact(x);
    error(k) = sqrt(h*sum((u -u_ref').^2));

    N = N*2;
end

%find the slope of the line using least squares fit
c = polyfit(log(h_list), log(error), 1)
log_err = polyval(c, log(h_list));
fprintf('Least square fit convergence rate = %f\n',c(1))




%plot the convergence plot
loglog(h_list,error,'o','LineWidth',2,'MarkerSize',10)
hold on
p = loglog(h_list,exp(log_err),'--k')
txt = ['slope = ',num2str(c(1))];
text(1.1*h_list(end),error(end),txt,'FontSize',15); %put text label on the figure
hold off
grid on;
xlabel('h'); ylabel('L2 error')
set(gca,'FontSize',18)


 
