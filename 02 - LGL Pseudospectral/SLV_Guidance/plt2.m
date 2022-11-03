%%

close all;

load Converged

sol2 = Optimal_Control_Problem();

%%
close all;

t2  = sol2{1};
x2  = sol2{2};
y2  = sol2{3};
vx2 = sol2{4};
vy2 = sol2{5};
u2  = sol2{6};

lamvx = lambda.eqnonlin(p.vxi)./p.w;
lamvy = lambda.eqnonlin(p.vyi)./p.w;

costates = atan2(-lamvy,-lamvx);

for i=1:length(t)
    xx1 = [x(i);y(i)];
    vv1 = [vx(i);vy(i)];
    gama(i) = acos(dot(xx1,vv1)/(norm(xx1)*norm(vv1)));
end

for i = 1:length(t2)
    xx2 = [x2(i);y2(i)];
    vv2 = [vx2(i);vy2(i)];
    gama2(i) = acos(dot(xx2,vv2)/(norm(xx2)*norm(vv2)));
end

%%
close all
figure(1);hold on;grid on
plot(t2,sqrt(x2.^2+y2.^2)/1e3-p.RE/1e3,'b.-');
plot(t,sqrt(x.^2+y.^2)/1e3-p.RE/1e3,'r--');
xlabel('t [sec]','FontSize',14,'FontWeight','Bold');
ylabel('h [km]','FontSize',14,'FontWeight','Bold');
title('Altitude','FontSize',14,'FontWeight','Bold');
legend('Indirect Method','Direct Method','FontSize',14,'FontWeight','Bold');

figure(2);hold on;grid on
plot(t2,sqrt(vx2.^2+vy2.^2),'b.-');
plot(t,sqrt(vx.^2+vy.^2),'r--');
xlabel('t [sec]','FontSize',14,'FontWeight','Bold');
ylabel('V [m/s]','FontSize',14,'FontWeight','Bold');
title('Velocity','FontSize',14,'FontWeight','Bold');
legend('Indirect Method','Direct Method','FontSize',14,'FontWeight','Bold');

figure(3);hold on;grid on
plot(t2,u2*180/pi,'b.-');
plot(t,u*180/pi,'r--');
plot(t,costates*180/pi-180,'g:');
xlabel('t [sec]','FontSize',14,'FontWeight','Bold');
ylabel('u [deg]','FontSize',14,'FontWeight','Bold');
title('Control Angle','FontSize',14,'FontWeight','Bold');
legend('Indirect Method','Direct Method','Direct Method Costate Estimates','FontSize',14,'FontWeight','Bold');

figure(4);hold on;grid on
plot(t2,(pi/2-gama2)*180/pi,'b.-');
plot(t,(pi/2-gama)*180/pi,'r--');
xlabel('t [sec]','FontSize',14,'FontWeight','Bold');
ylabel('\gamma [deg]','FontSize',14,'FontWeight','Bold');
title('Flight Parth Angle','FontSize',14,'FontWeight','Bold');
legend('Indirect Method','Direct Method','FontSize',14,'FontWeight','Bold');

figure(5);hold on;grid on
plot(t2,sqrt(vx2.^2+vy2.^2).*cos(gama2),'b.-');
plot(t,sqrt(vx.^2+vy.^2).*cos(gama'),'r--');
xlabel('t [sec]','FontSize',14,'FontWeight','Bold');
ylabel('Vr [m/s]','FontSize',14,'FontWeight','Bold');
title('Radial Velocity','FontSize',14,'FontWeight','Bold');
legend('Indirect Method','Direct Method','FontSize',14,'FontWeight','Bold');

figure(6);hold on;grid on
plot(t2,sqrt(vx2.^2+vy2.^2).*sin(gama2),'b--');
plot(t,sqrt(vx.^2+vy.^2).*sin(gama'),'r--');
xlabel('t [sec]','FontSize',14,'FontWeight','Bold');
ylabel('Vt [m/s]','FontSize',14,'FontWeight','Bold');
title('Tangential Velocity','FontSize',14,'FontWeight','Bold');
legend('Indirect Method','Direct Method','FontSize',14,'FontWeight','Bold');