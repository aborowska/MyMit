% cd ../
% cd ../
% addpath('include/');
% % addpath('include/NAIS');
% addpath('results/');
load('results\sv_VaR_IS_10000_CVtol_0.01.mat', 'theta_opt') 
 
M = size(theta_opt,1);

c = [theta_opt(1:M/2,1),theta_opt(1+M/2:M,1)];
phi = [theta_opt(1:M/2,2),theta_opt(1+M/2:M,2)];
sigma2 = [theta_opt(1:M/2,3),theta_opt(1+M/2:M,3)];
eta = [theta_opt(1:M/2,4),theta_opt(1+M/2:M,4)];
eps = [theta_opt(1:M/2,5),theta_opt(1+M/2:M,5)];


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
name = 'figures/presentation/sv_draws_param_both1.png';
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
name = 'figures/presentation/sv_draws_param_both2.png';
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
name = 'figures/presentation/sv_draws_param_both3.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
hold on
hist(eta,20)
h = findobj(gca, 'Type','patch');
set(h(1) ,'FaceColor', [0.7410 0.1 0.1], 'EdgeColor','w')
set(h(2) ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
% [f,xi] = ksdensity(c(:,1));
% plot(xi,f,'k');
% [f,xi] = ksdensity(c(:,2));
% plot(xi,f,'r');
% hold off
title('$$\eta_{\T+1}$$')
set(gca,'TickLabelInterpreter','latex')
name = 'figures/presentation/sv_draws_param_both4.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
hold on
hist(eps,20)
h = findobj(gca, 'Type','patch');
set(h(1) ,'FaceColor', [0.7410 0.1 0.1], 'EdgeColor','w')
set(h(2) ,'FaceColor', [0 0.4470 0.7410], 'EdgeColor','w')
% [f,xi] = ksdensity(c(:,1));
% plot(xi,f,'k');
% [f,xi] = ksdensity(c(:,2));
% plot(xi,f,'r');
% hold off
title('$$\eps_{\T+1}$$')
set(gca,'TickLabelInterpreter','latex')
name = 'figures/presentation/sv_draws_param_both5.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

% cd ../
% cd ../
% name = 'figures/presentation/sv_draws_param_both.png';
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(name,'-dpng','-r0')
