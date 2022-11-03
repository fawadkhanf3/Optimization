clear;close all;%clc

optimset fmincon;
opt = ans; %#ok<NOANS>
options = optimset(opt,'TolX',0.0001,'TolCon',0.0001,'MaxIter',5000);

N0  = 40;
sum = 0;
a = 0;
b = 1;
h = (b-a)/N0;

for i=1:N0+1
tt(i) = a+(i-1)*h;
end

N  = N0;
x1 = ones(N+1,1);
x2 = ones(N+1,1);
u  = zeros(N+1,1);
ud = zeros(N+1,1);
C2 = ones(N,1);
C3 = ones(N,1);

xx0 = [x1 ;x2; u ; ud ; C2; C3];

[xx,val,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@col_cost,xx0,[],[],[],[],[],[],@const_col2,options);

for i = 1:N+1
  
x1(i) = xx(i);
x2(i) = xx(N+1+i);
u(i)  = xx(2*(N+1)+i);
lamx(i) = lambda.eqnonlin(i);
lamz(i) = lambda.eqnonlin(N+1+i);

end

figure(1);hold on;
subplot(221);
plot(x1,x2,'r.-');grid on

subplot(222);
plot(tt,lamx,'r.-',tt,lamz,'b.-');grid on
legend('lamx','lamz');

subplot(223);
plot(tt,x1,'r.-',tt,x2,'b.-');grid on
legend('x','z');

subplot(224);
plot(tt,tan(u),'r.-',tt,lamz./lamx,'b.-');grid on
legend('tan(\gamma)','lamz/lamx');