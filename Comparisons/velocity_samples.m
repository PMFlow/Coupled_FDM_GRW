%% comparison of velocity samples:
%  for velocity fields computed by implicit and explicit FDM  

close all
I=201;
J=101;
a=0; b=20;
c=0; d=10;
dx=(b-a)/(I-1)
x=a:dx:b;
dy=(d-c)/(J-1)
y=c:dy:d;

figure; 

load ..\FDM_Flow\test_implicitFDM\dataGAUSS_100.mat
plot(x(1:199),Vx(:,50),x(1:199),Vy(:,50),'LineWidth',1.5); hold;
load ..\FDM_Flow\test_explicitFDM\dataGAUSS_100.mat
plot(x(1:199),Vx(:,50),'.',x(1:199),Vy(:,50),'.'); hold;
legend('$V_1(x,y=const)$','$V_2(x,y=const$)','Interpreter','latex','Location','best'); legend('boxoff');
xlabel('$x$','Interpreter','latex'); ylabel('Velocity samples computed by FDM');
title('implicit FDM (lines);  explcit FDM (dots)','FontWeight','normal','FontSize',10);
