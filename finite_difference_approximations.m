%%Finite difference approximations

%define a function
u = @(x) sin(x);

%define at which point we are approximating the derivative
x0 = 1;

%exact derivative for reference
dudx = @(x) cos(x)

%approximation formulations
Dp = @(h,x0) (u(x0+h) - u(x0))/h;
Dm = @(h,x0) (u(x0) - u(x0-h))/h;
D0 = @(h,x0) (u(x0+h) - u(x0-h))/(2*h);
D3 = @(h,x0) (2*u(x0+h) + 3*u(x0) - 6*u(x0-h) + u(x0 - 2*h))/(6*h);

%initial step size
h = 0.1

%repeat the experiment
for k=1:7
    error_Dp(k) = abs(Dp(h,x0) - dudx(x0));
    error_Dm(k) = abs(Dm(h,x0) - dudx(x0));
    error_D0(k) = abs(D0(h,x0) - dudx(x0));
    error_D3(k) = abs(D3(h,x0) - dudx(x0));

    %disp([h,x0,Dp(h,x0),Dm(h,x0),D0(h,x0),dudx(x0)])
    h_array(k) = h;
    h = h/3;
end

loglog(h_array,error_Dp,'-o',h_array,error_Dm,'--s',h_array,error_D0,'-.v',h_array,error_D3,':d','LineWidth',4)
grid on
legend('Dp','Dm','D0','D3')
xlabel('h'); ylabel('L2 error')
set(gca,'FontSize',18)

%calculate the ratio by which the error is decreasing
for k=2:5
    disp([h_array(k-1)/h_array(k),error_Dp(k-1)/error_Dp(k), error_Dm(k-1)/error_Dm(k), error_D0(k-1)/error_D0(k), error_D3(k-1)/error_D3(k)])
end