%% comparison of mean velocity and effective dispersion coeffiicents: 
%  for velocity fields computed by implicit and explicit FDM  

close all
D=0.01; 
%% mean velocity of center of mass 
figure;
load ensemble_coefficients_explicitFD_Nmod100R1e4
plot(t*dt,Mx,'r',t*dt,My,'b','LineWidth',1.5); hold all; ylim([-0.4 1.2]);
load linear_approximation
plot(t*dt,Mx,'.k',t*dt,My,'.k','MarkerSize',12);  
legend('$V_x(t)$','$V_y(t)$','Interpreter','latex','Location','best'); legend('boxoff');
xlabel('$t$','Interpreter','latex'); ylabel('Components of center of mass velocity'); set(gca, 'XTick', [0:2:10])
title('velocity by implicit FDM (lines);  linear approximation (dots)','FontWeight','normal','FontSize',10);
%% effective dispersion coefficients
figure
load ensemble_coefficients_explicitFD_Nmod100R1e4
plot(t*dt,Dx/D,'r',t*dt,Dy/D,'b','LineWidth',1.5); hold all;
load linear_approximation
plot(t*dt,Dx/D,'.k',t*dt,Dy/D,'.k','MarkerSize',12);
xlabel('$t$','Interpreter','latex'); ylabel('Dispersion coefficients'); set(gca, 'XTick', [0:2:10])
legend('$D_x(t) / D$','$D_y(t)/ D$','Interpreter','latex','Location','best'); legend('boxoff');
title('velocity by implicit FDM (lines);  linear approximation (dots)','FontWeight','normal','FontSize',10);


