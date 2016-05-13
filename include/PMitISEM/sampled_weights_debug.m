ind_fin = isfinite(lnk);
sum(ind_fin)
theta_fin = theta(ind_fin,:);
lnk_fin = lnk(ind_fin);
lnd_fin = lnd(ind_fin);

reshape(mit_init.Sigma,11,11)


w_fin = lnk_fin - lnd_fin;

[w_fin_sort, ind_sort] = sort(w_fin); 
theta_fin_sort = theta_fin(ind_sort,:);
lnk_fin_sort = lnk_fin(ind_sort);
lnd_fin_sort = lnd_fin(ind_sort);


y_fin_sort = predict_arch(theta_fin_sort(:,1), y_T, S, H, theta_fin_sort(:,2:end));
PL_fin_sort = (fn_PL(y_fin_sort));

figure(100)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1)    
    hold on
    plot(w_fin,'Color',[0.5 0.5 0.5])
    plot(w_fin_sort, 'b')
    title('w fin and fin sort')
    hold off
subplot(2,2,2)
    hist(w_fin_sort,20)
    title('histogram w fin sort')   
subplot(2,2,3)
    plot(PL_fin_sort)
    title('PL fin sort')
subplot(2,2,4)
    scatter(w_fin_sort,PL_fin_sort)
    xlabel('w fin sort') 
    ylabel('PL fin sort')
    
name = ['figures/PMitISEM/',model,'_H', num2str(H),'_sampled_weights.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')
        

figure(200)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
subplot(1,3,1)
plot(theta_fin_sort(1:10,2:end)')
subplot(1,3,2)
plot(theta_fin_sort(4000-10:4000,2:end)')
subplot(1,3,3)
plot(theta_fin_sort(end-10:end,2:end)')
 

figure(300)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
    [PL_sort, ind_PL] = sort(PL_fin_sort);
    plot(PL_sort)
    title('PL sort')
subplot(1,2,2)
    plot(w_fin_sort(ind_PL))
    title('w corresp. to PL sort')
name = ['figures/PMitISEM/',model,'_H', num2str(H),'_PL_sort.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')

theta_high_weight = theta_fin_sort(end-10:end,:);
y_high_weight = predict_arch(theta_high_weight(:,1), y_T, S, H, theta_high_weight(:,2:end));
PL_high_weight = fn_PL(y_high_weight);
w_high_weight = kernel(theta_high_weight) - dmvgt(theta_high_weight, mit_init, true, GamMat);

lnk_high_weight = kernel(theta_high_weight);
lnd_high_weight = dmvgt(theta_high_weight, mit_init, true, GamMat);




theta_median = theta_fin_sort(4000-10:4000,:);
y_median = predict_arch(theta_median(:,1), y_T, S, H, theta_median(:,2:end));
PL_median = fn_PL(y_median);
w_median = kernel(theta_median) - dmvgt(theta_median, mit_init, true, GamMat);

lnk_median = kernel(theta_median);
lnd_median = dmvgt(theta_median, mit_init, true, GamMat);

figure(400)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
hist([lnk_high_weight,lnk_median],22)
legend('lnk high weight','lnk median')
subplot(1,2,2)
hist([lnd_high_weight,lnd_median],22)
legend('lnd high weight','lnd median')
name = ['figures/PMitISEM/',model,'_H', num2str(H),'_lnd_lnk.png'];
set(gcf,'PaperPositionMode','auto');
print(name,'-dpng','-r0')



ind_sort_trunc = ind_sort(round(0.1*length(ind_sort)):round(0.9*length(ind_sort)));
w_fin_sort_trunc = w_fin(ind_sort_trunc); 
theta_fin_sort_trunc = theta_fin(ind_sort_trunc,:);

w_fin_sort_trunc = w_fin_sort_trunc-max(w_fin_sort_trunc);
w_fin_sort_trunc=exp(w_fin_sort_trunc);
