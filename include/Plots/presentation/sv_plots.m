cd ../
cd ../
addpath('include/');
% addpath('include/NAIS');
addpath('results/');
load('sv_mitisem.mat', 'theta1', 'theta_hl') 


c = [theta1(:,1), theta_hl(:,1)];
phi = [theta1(:,2), theta_hl(:,2)];
sigma2 = [theta1(:,3), theta_hl(:,3)];

figure(100)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');

subplot(3,1,1)
hold on
hist(c,20)
h = findobj(gca, 'Type','patch');
set(h(1) ,'FaceColor', [0.7410 0.1 0.1], 'EdgeColor','w')
set(h(2) ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
% [f,xi] = ksdensity(c(:,1));
% plot(xi,f,'k');
% [f,xi] = ksdensity(c(:,2));
% plot(xi,f,'r');
% hold off
title('$$c$$')
set(gca,'TickLabelInterpreter','latex')

subplot(3,1,2)
hold on
hist(phi,20)
h = findobj(gca, 'Type','patch');
set(h(1) ,'FaceColor', [0.7410 0.1 0.1], 'EdgeColor','w')
set(h(2) ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
% [f,xi] = ksdensity(c(:,1));
% plot(xi,f,'k');
% [f,xi] = ksdensity(c(:,2));
% plot(xi,f,'r');
% hold off
title('$$\phi$$')
set(gca,'TickLabelInterpreter','latex')

subplot(3,1,3)
hold on
hist(sigma2,20)
h = findobj(gca, 'Type','patch');
set(h(1) ,'FaceColor', [0.7410 0.1 0.1], 'EdgeColor','w')
set(h(2) ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
% [f,xi] = ksdensity(c(:,1));
% plot(xi,f,'k');
% [f,xi] = ksdensity(c(:,2));
% plot(xi,f,'r');
% hold off
title('$$\sigma^2_{\eta}$$')
set(gca,'TickLabelInterpreter','latex')


cd ../
cd ../
name = 'figures/presentation/sv_draws_param_both.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')
