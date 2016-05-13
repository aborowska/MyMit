
figure(150)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1)
kernel = @(xx) posterior_arch(xx, data, S, true);
logpdf = kernel((0:0.01:0.99)');
plot((0:0.01:0.99)',logpdf)
xlabel('\alpha')
title('Log likelihood')

subplot(2,2,2)
lognorm = - 0.5*(log(2*pi) + (-10:0.1:10).^2);
plot((-10:0.1:10),lognorm)
title('Log of the standard normal density')

subplot(2,2,3)
alpha_lnk = kernel(draw_hl(:,1));
[sort_alpha_lnk, ind] = sort(alpha_lnk);
alpha_sort = draw_hl(ind, 1);
plot(alpha_sort, sort_alpha_lnk, 'LineStyle','None','Marker','.')
xlabel('high loss \alpha draws')
title('Log likelihood of high loss \alpha draws')


subplot(2,2,4)
logeps = sum(- 0.5*(log(2*pi) + (draw_hl(:,2:11)).^2),2);
plot(logeps)
title('Log density of  10-dim. high loss \epsilon draws')

hold on
plot(logeps/10,'r')
plot(alpha_lnk/T)

hold off

figure(129)

