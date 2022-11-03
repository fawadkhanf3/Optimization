function J = obj_fun(x,G)

wp.M  = x(1);
wp.A  = x(2);
wp.wB = x(3);

Wp = tf([1/wp.M wp.wB],[1 wp.wB*wp.A]);

wt.M  = x(4);
wt.A  = x(5);
wt.wB = x(6);

Wt = tf([1/wt.M wt.wB],[1 wt.wB*wt.A]);

Wu = tf(x(7));

try
    [K,CL,GAM,INFO] = mixsyn(G,Wp,Wu,Wt);
    
    L = K*G;
    S = 1/(1+L);
    T = 1-S;
    
    J = norm([Wp*S;Wt*T;Wu*K*S],'inf');
    
    
    
    
catch
    J = 1e9;
end

fprintf('\n J = %g',J);

end