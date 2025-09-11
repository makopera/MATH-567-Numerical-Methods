clear
%forcing function
f = @(x) exp(x);

sigma = 0;
b = 2;
xa = 0;
xb = 1;

L = xb - xa;  %length of the domain

%exact solution
u_exact = @(x) exp(x) + (sigma-exp(xa))*(x - xb) + b - exp(xb);

%% We will solve d^2(u)/dx^2 = f(x) on domain x = [x_a,x_b] 
% with BC: 
%   u'(x_a) = sigma
%   u(x_b)  = b
% using a central finite difference method

neumann_method= 1; %1 - first order one sided difference

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
    case 2 
        A(1,1) = -h;
        A(1,2) = h;
end

%modify last row to represent BC
A(M,M-1) = 0;
A(M,M) = h^2;

A = A/h^2;

%NOTE: for large matrices it is a better idea to use sparse matrix
%representation

%create RHS vector
g = (f(x))';

%implement BC
switch neumann_method
    case 1
        g(1) = sigma;
    case 2
        g(1) = sigma + h/2*f(x(1));
end



g(M) = b;

%solve linear system
u = A\g;

%plot solution
xx= linspace(xa,xb);
plot(xx,u_exact(xx),'--k','LineWidth',2)
hold on;
plot(x,u,'o','LineWidth',2,'MarkerSize',10);
%plot(xa,a,'ob','LineWidth',2,'MarkerSize',10);
%plot(xb,b,'ob','LineWidth',2,'MarkerSize',10);

hold off
xlabel('x'); ylabel('u(x)');
grid on
set(gca,'FontSize',18)

%% Convergence study

neumann_method= 2; %1 - first order one sided difference
N = 5;
n_tests = 5;
error = zeros(1,n_tests);
h_list = zeros(1,n_tests); 

for k=1:n_tests
    
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
        case 2 
            A(1,1) = -h;
            A(1,2) = h;
    end

    %modify last row to represent BC
    A(M,M-1) = 0;
    A(M,M) = h^2;

    A = A/h^2;

    %NOTE: for large matrices it is a better idea to use sparse matrix
    %representation

    %create RHS vector
    g = (f(x))';

    %implement BC
    switch neumann_method
        case 1
            g(1) = sigma;
        case 2
            g(1) = sigma + h/2*f(x(1));
    end
    g(M) = b;

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

%% Compare all the methods
 
for neumann_method = 1:2

    N = 5;
    n_tests = 5;
    error = zeros(1,n_tests);
    h_list = zeros(1,n_tests);

    for k=1:n_tests

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
            case 2
                A(1,1) = -h;
                A(1,2) = h;
        end

        %modify last row to represent BC
        A(M,M-1) = 0;
        A(M,M) = h^2;

        A = A/h^2;

        %NOTE: for large matrices it is a better idea to use sparse matrix
        %representation

        %create RHS vector
        g = (f(x))';

        %implement BC
        switch neumann_method
            case 1
                g(1) = sigma;
            case 2
                g(1) = sigma + h/2*f(x(1));
        end
        g(M) = b;

        %solve linear system
        u = A\g;

        %calculate error
        h_list(k) = h;
        u_ref = u_exact(x);
        error(k) = sqrt(h*sum((u -u_ref').^2));

        N = N*2;
    end

    %find the slope of the line using least squares fit
    c = polyfit(log(h_list), log(error), 1);
    

    error_methods(neumann_method,:) = error;
    c_methods(neumann_method,:) = c;

    



end

%plot the convergence plot
    loglog(h_list,error_methods(1,:),'o','LineWidth',2,'MarkerSize',10)
    hold on
    log_err = polyval(c_methods(1,:), log(h_list));
    fprintf('Least square fit convergence rate = %f\n',c_methods(1,1))
    p = loglog(h_list,exp(log_err),'--k')
    txt = ['slope = ',num2str(c_methods(1,1))];
    text(1.1*h_list(end),error_methods(1,end),txt,'FontSize',15); %put text label on the figure

    loglog(h_list,error_methods(2,:),'o','LineWidth',2,'MarkerSize',10)

    log_err = polyval(c_methods(2,:), log(h_list));
    fprintf('Least square fit convergence rate = %f\n',c_methods(2,1))
    p = loglog(h_list,exp(log_err),'--k')
    txt = ['slope = ',num2str(c_methods(2,1))];
    text(1.1*h_list(end),error_methods(2,end),txt,'FontSize',15); %put text label on the figure


    hold off
    grid on;
    xlabel('h'); ylabel('L2 error')
    set(gca,'FontSize',18)