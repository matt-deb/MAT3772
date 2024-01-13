clear;
clc;

%order of the system
M=4;

%setting step size, start and end times, and number of mesh points
h=0.01;
t0=0;
tN=1;
N=(tN-t0)/h;

%creating mesh
T = linspace(t0,tN,N+1);

%defining f(t,y)
F = @(t,y) [y(2,:); y(3,:); y(4,:); y(1,:)+2*t];

%initializing Y and setting initial conditions
Y = zeros(M,N+1);
Y(:,1) = [1;-2;1;2];

%calculating y_1 using predictor-corrector method as in the for loop
y_p = Y(:,1) + h*F(T(1),Y(:,1));
Y(:,2) = Y(:,1) + h*F(T(1),Y(:,1));

for n=2:N
    %predictor: Euler method
    y_p = Y(:,n) + h*F(T(n),Y(:,n));

    %corrector: Simpson's method
    Y(:,n+1) = Y(:,n-1) + (h/3)*(F(T(n-1),Y(:,n-1))+4*F(T(n),Y(:,n))+F(T(n+1),y_p));
end

%defining exact solution
y_exc = [exp(T) - sin(T) - 2*T; exp(T) - cos(T) - 2; exp(T) + sin(T); exp(T) + cos(T)];

%plotting
plot(T,y_exc(1,:),'r-')
hold on
plot(T,Y(1,:),'b o','MarkerFaceColor','b','MarkerSize',1)
hold off

%finding relative error and displaying result
diff = y_exc(:,N+1) - Y(:,N+1);
rel_err = sum(diff.^2,1)/sum(y_exc(:,N+1).^2,1);
fprintf("The relative error at t_N = 1 is %d.", rel_err)
