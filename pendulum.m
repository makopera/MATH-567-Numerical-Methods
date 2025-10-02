%% We are solving a pendulum equation d^2 theta / dt^2 = -g/L*sin(theta)

clear
clf

g = 10;
L = 10;


%boudnary conditions - assume Dirichlet at the beginning and the end
a = pi;
b = 0.7;
T = 2*pi;


%% Start with linearized equation 

N = 20; %number of internal gridpoints -  here time points
h = T/(N+1); %interval length
M=N+2; %account for all points including boundary

%create grid 
t = linspace(0,T,M); %this mesh includes endpoints

%build matrix A
A = zeros(M,M); %initialize zero matrix N+2 x N+2 
A = diag(-2*ones(M,1)+g/L*h^2)    + ...    % main diagonal
    diag(ones(M-1,1), -1) + ...    % upper diagonal
    diag(ones(M-1,1),  1);         % lower diagonal

%modify first row to represent BC
A(1,1) = h^2; %Dirichlet
A(1,2) = 0;
%modify last row to represent BC
A(M,M-1) = 0;
A(M,M) = h^2; %Dirichlet

A = A/h^2

%NOTE: for large matrices it is a better idea to use sparse matrix
%representation

%create RHS vector
g = zeros(M,1);

%implement BC
g(1) = a;
g(M) = b;

%solve linear system
theta = A\g;

%plot solution
%plot(t,theta,'--k','LineWidth',2,'MarkerSize',10);

%hold on
%xlabel('t'); ylabel('\theta(t)');
%grid on
%set(gca,'FontSize',18)

%% Non-linear equation

N = 20; %number of internal gridpoints -  here time points
h = T/(N+1); %interval length
M=N+2; %account for all points including boundary

%create grid 
t = linspace(0,T,M); %this mesh includes endpoints

%initial guess
%theta =  0*ones(1,M); %0.7*cos(t) + 0.5*sin(t); %a*ones(1,M);
theta = a*ones(1,M) + sin(t/2);
theta(1) = a;
theta(M) = b;
theta = theta';

[G,J] = pendulum_matrix(theta,N,h);

plot(t,theta,":",'LineWidth',2)
hold on

for k=1:9
    [G,J] = pendulum_matrix(theta,N,h);
    delta = J\(-G);
    theta(2:end-1) = theta(2:end-1) + delta;
    plot(t,theta,":",'LineWidth',2)
    fprintf("k=%d, norm(delta) = %e\n",k,max(delta))
end

plot(t,theta,"-.b",'LineWidth',2)
hold off
legend('linearised','k=1','k=2','k=3','k=4','k=5','k=6','k=7','k=8','k=9','k=10')

function [G,J] = pendulum_matrix(theta,N,h)

    %build vector G(theta)
    G = zeros(N,1); %initialize zero vector
    for i=2:N+1
        G(i-1) = 1/h^2*(theta(i-1) - 2*theta(i) + theta(i+1))+sin(theta(i));
    end

    %build Jacobian matrix J_ij = \frac{\partial G_i(\theta)}{\partial theta_j}
    
    J = diag(-2*ones(N,1)+h^2*cos(theta(2:end-1)))    + ...    % main diagonal
    diag(ones(N-1,1), -1) + ...    % upper diagonal
    diag(ones(N-1,1),  1);         % lower diagonal
    J = J/h^2;
    
    % for i=1:m
    %     for j=1:m
    %         if((j==i-1)||(j=i+1))
    %             J(i,j) = 1/h^2
    %         elseif j==1
    %             J(i,j) = 
    %         end
    %     end
    % end
end