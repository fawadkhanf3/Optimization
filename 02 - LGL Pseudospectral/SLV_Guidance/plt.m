%% 

close all;

load Converged

for i=1:length(t)
    xx1 = [x(i);y(i)];
    vv1 = [vx(i);vy(i)];
    gama(i) = acos(dot(xx1,vv1)/(norm(xx1)*norm(vv1)));
end

figure(1);hold on;grid on
plot(t,sqrt(x.^2+y.^2)/1e3-p.RE,'r.-');
xlabel('t [sec]');ylabel('h [km]');
title('Altitude');

figure(2);hold on;grid on
plot(t,sqrt(vx.^2+vy.^2),'r.-');
xlabel('t [sec]');ylabel('V [m/s]');
title('Velocity');

figure(3);hold on;grid on
plot(t,u*180/pi,'r.-');
xlabel('t [sec]');ylabel('u [deg]');
title('Control Angle');

figure(4);hold on;grid on
plot(t,(pi/2-gama)*180/pi,'r.-');
xlabel('t [sec]');ylabel('\gamma [deg]');
title('Flight Parth Angle w.r.t Local Horizontal');

figure(5);hold on;grid on
plot(t,sqrt(vx.^2+vy.^2).*cos(gama'),'r.-');
xlabel('t [sec]');ylabel('Vr [m/s]');
title('Radial Velocity');

figure(6);hold on;grid on
plot(t,sqrt(vx.^2+vy.^2).*sin(gama'),'r.-');
xlabel('t [sec]');ylabel('Vt [m/s]');
title('Tangential Velocity');

lamx = lambda.eqnonlin(p.xi);
lamy = lambda.eqnonlin(p.yi);
lamvx = lambda.eqnonlin(p.vxi);
lamvy = lambda.eqnonlin(p.vyi);

u2 = atan2(-lamvy,-lamvx);

figure(7);hold on;grid on
plot(t,u2*180/pi,'r.-');
xlabel('t [sec]');ylabel('u2 [deg]');
title('Control Angle');