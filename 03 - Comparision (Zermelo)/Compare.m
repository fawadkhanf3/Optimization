%% Comparison


clear;close all;clc

%%

yCollocation    = Zermelo_Collocation();
yPseudospectral = Zermelo_Pseudospectral();

%%
figure(1);hold on;grid on;
plot(yCollocation.x,yCollocation.z,'r.-');
plot(yPseudospectral.x,yPseudospectral.z,'b--');
xlabel('x');ylabel('z');
legend('Collocation','Pseudospectral','Location','Best');

figure(2);hold on;grid on;
plot(yCollocation.t,yCollocation.x,'r.-');
plot(yPseudospectral.t,yPseudospectral.x,'b--');

plot(yCollocation.t,yCollocation.z,'g.-');
plot(yPseudospectral.t,yPseudospectral.z,'m--');

xlabel('t');ylabel('States');
legend('x (Collocation)','x (Pseudospectral)',...
    'z (Collocation)','z (Pseudospectral)',...
    'Location','Best');

figure(3);hold on;grid on;
plot(yCollocation.t,tan(yCollocation.u),'r.-');
plot(yCollocation.t(1:end-1),yCollocation.lamz./yCollocation.lamx,'b--');

plot(yPseudospectral.t,tan(yPseudospectral.u),'g.-');
plot(yPseudospectral.t,yPseudospectral.lamz./yPseudospectral.lamx,'m--');

xlabel('t');ylabel('Control');
legend('tan(\gamma) (Collocation)','\lambda_{z}/\lambda_{x} (Collocation)',...
    'tan(\gamma) (Pseudospectral)','\lambda_{z}/\lambda_{x} (Pseudospectral)',...
    'Location','Best');

figure(4);hold on;grid on;
plot(yCollocation.t(1:end-1),yCollocation.lamx,'r.-');
plot(yPseudospectral.t,yPseudospectral.lamx,'b--');

plot(yCollocation.t(1:end-1),yCollocation.lamz,'g.-');
plot(yPseudospectral.t,yPseudospectral.lamz,'m--');

xlabel('t');ylabel('CoStates');
legend('\lambda_{x} (Collocation)','\lambda_{x} (Pseudospectral)',...
    '\lambda_{z} (Collocation)','\lambda_{z} (Pseudospectral)',...
    'Location','Best');