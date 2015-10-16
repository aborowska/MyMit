clear all

cd ../
cd ../
addpath('include/');
% addpath('include/NAIS');
addpath('results/');
load('svt_mitisem.mat', 'theta1', 'theta_hl') 

c = [theta1(:,1), theta_hl(:,1)];
phi = [theta1(:,2), theta_hl(:,2)];
sigma2 = [theta1(:,3), theta_hl(:,3)];
nu = [theta1(:,4), theta_hl(:,4)];


figure(1)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
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
name = 'figures/presentation/svt_draws_param_both1.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
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
name = 'figures/presentation/svt_draws_param_both2.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
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
name = 'figures/presentation/svt_draws_param_both3.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
hold on
hist(nu,20)
h = findobj(gca, 'Type','patch');
set(h(1) ,'FaceColor', [0.7410 0.1 0.1], 'EdgeColor','w')
set(h(2) ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
% [f,xi] = ksdensity(c(:,1));
% plot(xi,f,'k');
% [f,xi] = ksdensity(c(:,2));
% plot(xi,f,'r');
% hold off
title('$$\nu$$')
set(gca,'TickLabelInterpreter','latex')
name = 'figures/presentation/svt_draws_param_both4.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

