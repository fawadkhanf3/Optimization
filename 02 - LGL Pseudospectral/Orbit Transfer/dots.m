function dx = dots(t,x,mm,p)

r = x(1);
theta = x(2);
u = x(3);
v = x(4);

E = interp1(mm(:,1),mm(:,2),t);

dr     = u;
dtheta = v./r;
du     = (v.^2)./r - 1./(r.^2) + p.A.*sin(E);
dv     = -(u.*v)./r + p.A.*cos(E);

dx = [dr;dtheta;du;dv];

end