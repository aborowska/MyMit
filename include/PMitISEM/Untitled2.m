%% loglik
loglik = zeros(M,1);
alpha = draw_hl(:,1);
omega = S*(1-alpha);
for ii = 1:M
%     h = S;
%     loglik(ii,1) = -0.5*(log(2*pi) + log(h) + (data(1)^2)/h); 
     for jj = 2:T
        h = omega(ii,1) + alpha(ii,1)*data(jj-1)^2;
        loglik(ii,1) = loglik(ii,1) -0.5*(log(2*pi) + log(h) + (data(jj)^2)/h); 
    end
    if ((alpha(ii,1) >= 0) & (alpha(ii,1) < 1))
        loglik(ii,1) = loglik(ii,1) + 1;
    else
        loglik(ii,1) = loglik(ii,1) - Inf;
    end
end

kernel = @(xx) posterior_arch(xx, data, S, true);
loglik_kernel = kernel(alpha);


loglik_eps = -0.5*(H*log(2*pi) + sum(draw_hl(:,2:H+1).^2,2));


%% get draws and log kernel evaluation

xxx = -0.4:0.01:1;
xxx = xxx';
lnk_xxx = kernel(xxx);
lnd_xxx = dmvgt(xxx, mit_one, true, GamMat);
w_xxx=  fn_ISwgts(lnk_xxx, lnd_xxx, norm);

lnd_xxx2 = dmvgt(xxx, mit_new, true, GamMat);
w_xxx2 =  fn_ISwgts(lnk_xxx, lnd_xxx2, norm);

figure(111)
subplot(2,2,1)
plot(xxx,lnk_xxx)

subplot(2,2,2)
hold on
plot(xxx,lnd_xxx)
plot(xxx,lnd_xxx2,'r')
hold off

subplot(2,2,3)
hold on
plot(xxx,w_xxx)
plot(xxx,w_xxx2,'r')
hold off  

subplot(2,2,4)
%     hold all
%     plot(xxx,lnk_xxx)
%     plot(xxx,lnd_xxx)
%     plot(xxx,w_xxx)
%     hold off
k_xxx = exp(lnk_xxx-max(lnk_xxx));
plot(xxx,k_xxx)

d_xxx = exp(lnd_xxx-max(lnd_xxx));
d_xxx2 = exp(lnd_xxx2-max(lnd_xxx2));

figure(2222)
hold on 
plot(xxx,k_xxx)
plot(xxx,d_xxx,'g')
plot(xxx,d_xxx2,'r')
hold off

    
 %% draw_hl hist   
    
figure(456)
set(gcf,'units','normalized','outerposition',[0 0 0.3 1]);
subplot(3,1,1)
hist(draw_hl(:,1))
title('draw_h_l: \alpha')
subplot(3,1,2)
hist(draw_hl(:,2))
title('draw_h_l: \epsilon_1')
subplot(3,1,3)
hist(draw_hl(:,3))
title('draw_h_l: \epsilon_2')
suptitle(['ARCH H=',num2str(H),' draw_h_l'])

name = ['figures/PMitISEM/arch_H',num2str(H),'_draw_hl_hist.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')


figure(4567)
set(gcf,'units','normalized','outerposition',[0 0 0.3 1]);
subplot(3,1,1)
hist(theta(:,1))
title('theta: \alpha')
subplot(3,1,2)
hist(theta(:,2))
title('theta: \epsilon_1')
subplot(3,1,3)
hist(theta(:,3))
title('theta: \epsilon_2')
suptitle(['ARCH H=',num2str(H),' theta from mit init'])

name = ['figures/PMitISEM/arch_H',num2str(H),'_theta_init_hist.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')



figure(45678)
set(gcf,'units','normalized','outerposition',[0 0 0.3 1]);
subplot(3,1,1)
hist([draw_hl(:,1),theta(:,1)])
title('draw_h_l, theta: \alpha')
subplot(3,1,2)
hist([draw_hl(:,2),theta(:,2)])
title('draw_h_l, theta: \epsilon_1')
subplot(3,1,3)
hist([draw_hl(:,3),theta(:,3)])
title('draw_h_l, theta: \epsilon_2')
suptitle(['ARCH H=',num2str(H),' draw_h_l and theta from mit init'])

name = ['figures/PMitISEM/arch_H',num2str(H),'_draw_hl_theta_init_hist.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')


loglik_eps = -0.5*(H*log(2*pi) + sum(draw_hl(:,2:H+1).^2,2)); 

figure(12)
set(gcf,'units','normalized','outerposition',[0 0 0.3 1]);
subplot(3,1,1)
plotyy(1:M,lnk_hl,1:M,loglik_eps)
title('lnk_h_l (on posterior) and loglik of \epsilon(s) (sorted wrt PL)')
subplot(3,1,2)
plot(lnd_hl)
title('lnd_h_l (on mit1 for posterior)')
subplot(3,1,3)
plot(w_hl)
title('w_h_l')
suptitle(['ARCH H=',num2str(H),' w_h_l and components (only of \alpha, for \epsilon is 1)'])

name = ['figures/PMitISEM/arch_H',num2str(H),'_w_draw_hl.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')



loglik_eps_init = -0.5*(H*log(2*pi) + sum(theta(:,2:H+1).^2,2)); 

figure(123)
set(gcf,'units','normalized','outerposition',[0 0 0.3 1]);
subplot(3,1,1)
plotyy(1:M,lnk,1:M,loglik_eps_init )
title('lnk (on posterior + target for \epsilon(s)), and loglik for  \epsilon(s))')
subplot(3,1,2)
plot(lnd)
title('lnd (on mix init for \theta)')
subplot(3,1,3)
plot(w)
title('w')
suptitle(['ARCH H=',num2str(H),' w for mit init and components'])

name = ['figures/PMitISEM/arch_H',num2str(H),'_w_theta_init.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')



figure(1234)
set(gcf,'units','normalized','outerposition',[0 0 0.3 1]);
subplot(3,1,1)
plot(lnk)
title('lnk (on posterior + target for \epsilon(s))')
subplot(3,1,2)
plot(lnd)
title('lnd (on mix init for \theta)')
subplot(3,1,3)
plot(lnk-lnd)
title('lnk-lnd')
suptitle(['ARCH H=',num2str(H),' w for mit init and components'])

name = ['figures/PMitISEM/arch_H',num2str(H),'_w_theta_init2.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')



%% max lnw?
lnw = lnk-lnd;
[lnw_sort,ind_w] = sort(lnw);
loglik_eps_init_sort = loglik_eps_init(ind_w);

y_init = predict_arch(theta(:,1), y_T, S, H, theta(:,2:H+1));  
PL_init = fn_PL(y_init);
PL_init_sort = PL_init(ind_w);

figure(1234)
set(gcf,'units','normalized','outerposition',[0 0 0.3 1]);
subplot(3,1,1)
plot(lnw,'color',[0.5 0.5 0.5])
hold on
[hAx,hLine1,hLine2] = plotyy(1:M,lnw_sort,1:M,loglik_eps_init_sort)
set(hLine1,'LineWidth',2)
hold off
title('lnw = lnk-lnd, lnw sorted, loglik for \epsilon(s) sorted')
subplot(3,1,2)
plotyy(1:M,theta(ind_w,1),1:M,[theta(ind_w,2),theta(ind_w,3)])
title('\theta sort')
legend('\alpha','\epsilon_1','\epsilon_2','Location','SouthEast')
subplot(3,1,3)
plot(PL_init_sort)
title('PL sort')
suptitle(['ARCH H=',num2str(H),' lnw sort'])

name = ['figures/PMitISEM/arch_H',num2str(H),'_lnw_theta_PL_init_sort.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')



figure(444)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.75]);

plotyy(1:M,PL_init_sort, 1:M, lnw_sort)
hold on
plot(1:M,ones(M,1)*mean(VaR_prelim),'r','LineWidth',2)
hold off
title(['ARCH H=',num2str(H),' PL init sort, lnw sort and VaR prelim'])
name = ['figures/PMitISEM/arch_H',num2str(H),'_PL_init_sort.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')

 