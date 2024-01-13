clear;
clc;

%setting values of h and start and end times, and initialzing array to hold
%errors
H = [1/300,1/3000];
t0 = 0;
tN = 2/3;
err = zeros(2,1);

%defining rhs of IVP and exact solution
F = @(t,y) t*y/(1-t^2);
y_exc = @(t) -(sqrt(1-t.^2)).^(-1);

for i=1:2
    %choosing value of h
    h = H(i);

    %calculating number of mesh points
    N = int32((tN-t0)/h);

    %creating mesh
    T = linspace(t0,tN,N+1);

    %creating array to hold calculations of each y and setting initial
    %condition
    Y = zeros(N+1,1);
    Y(1) = -1;

    %calculating y_1 using the modified Euler method
    f0 = F(T(1),Y(1));
    Y(2) = Y(1) + h*F(T(1)+h/2,Y(1)+(h/2)*f0);

    %calculating other y's using the midpoint rule
    for n=1:N-1
        fn = F(T(n),Y(n));
        Y(n+2) = Y(n) + 2*h*F(T(n+1),Y(n+1));  
    end

    %plotting numerical and exact solutions
    plot(T,y_exc(T),'r-')
    hold on
    plot(T,Y,'b o','MarkerFaceColor','b','MarkerSize',1)
    hold off

    %calculating error
    err(i) = norm(Y(end)-y_exc(T(end)));
end

%calculating order of convergence and printing result
p = log(err(1)/err(2))/log(H(1)/H(2));
fprintf("The order of convergence is %f.", p)