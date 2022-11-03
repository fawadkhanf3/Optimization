function [c,d] = const_col2(xx)

N = 40;
a = 0;
b = 1;
h = (b-a)/N;

for i=1:N+1
    x1(i) = xx(i);
    x2(i) = xx(N+1+i);
    u(i)  = xx(2*(N+1)+i);
    ud(i) = xx(3*(N+1)+i);
end

for i=1:N
    C2(i) = xx(4*(N+1)+i);
    C3(i) = xx(4*(N+1)+N+i);
    C1(i) = ud(i);
    C0(i) = u(i);
end

for i=1:N+1
    f1(i) = 0.3*x2(i)+cos(u(i));
    f2(i) = sin(u(i));
end

for i=1:N
    x1c(i) = (x1(i)+x1(i+1))/2+h*(f1(i)-f1(i+1))/8;
    x2c(i) = (x2(i)+x2(i+1))/2+h*(f2(i)-f2(i+1))/8;
    
    uc(i)  = C0(i)+C1(i)/2+C2(i)/4+C3(i)/8;
    
    fc1(i) = 0.3*x2c(i)+cos(uc(i));
    fc2(i) = sin(uc(i));
   
    x1d(i) = -3*(x1(i)-x1(i+1))/2/h -(f1(i)+f1(i+1))/4;
    x2d(i) = -3*(x2(i)-x2(i+1))/2/h -(f2(i)+f2(i+1))/4;
    
    d1(i)  = fc1(i)-x1d(i);
    d2(i)  = fc2(i)-x2d(i);
end

for i=1:N
    d(i)   = d1(i);
    d(N+i) = d2(i);
end

for i=1:N
    d(2*N+i) = C0(i)+C1(i)+C2(i)+C3(i)-u(i+1);
    d(3*N+i) = C1(i)+2*C2(i)+3*C3(i)-ud(i+1);
end

for i=1:N-1
    d(4*N+i) = 2*C2(i)+6*C3(i)-2*C2(i+1);
end

d(5*N)   = x1(1)-0;
d(5*N+1) = x2(1)-0;

c = [u-pi -pi-u];