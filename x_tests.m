ACF = zeros(M,1);

for ii=1:M
    pom = autocorr(x(ii,:),2);
    ACF(ii,1) = pom(2);
end
    
    f_pl = @(aa) 100*(exp(aa/100) - 1); 

    ind_y = (imag(y_opt_h1)==0); 
    y_opt_h1 = y_opt_h1(ind_y); 
    PL_opt_h1 = f_pl(sum(y_opt_h1,2));
    [PL_opt_h1, ind] = sort(PL_opt_h1); 
    
    pos =  max(find(PL_opt_h1<= VaR_IS(sim)));
 
    theta_opt_sort = theta_opt(ind_y,:);
    theta_opt_sort = theta_opt_sort(ind,:);
    
    lnk_sort = lnk(ind_y,:);
    lnk_sort = lnk_sort(ind,:);
 
    x_sort = x(ind_y,:);
    x_sort = x_sort(ind,:);
    
    lng_y_sort = lng_y(ind_y,:);
    lng_y_sort = lng_y_sort(ind,:);
    
    lnw_x_sort = lnw_x(ind_y,:);
    lnw_x_sort = lnw_x_sort(ind,:);
    
    eps_bar_sort = eps_bar(ind_y,:);
    eps_bar_sort = eps_bar_sort(ind,:);
    
    eps_sim_sort = eps_sim(ind_y,:);
    eps_sim_sort = eps_sim_sort(ind,:);
    
    C_sim_sort = C_sim(ind_y,:);
    C_sim_sort = C_sim_sort(ind,:);
    
    lnp_T_sort = lnp_T(ind_y,:);
    lnp_T_sort = lnp_T_sort(ind,:);
    
    y_star_sort = y_star(ind_y,:);
    y_star_sort = y_star_sort(ind,:);
 
figure(333)
hold all
for ii = 1:5
    q = y_star_sort((abs(y_star_sort(:,ii))<100),ii);
    plot(q)
end
hold off
  
figure(22)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');   
for ii = 1:4
subplot(2,2,ii)
hold on
q = lnp_T_sort((imag(lnp_T_sort(:,end-4+ii))==0),end-4+ii);
plot(q)
line([pos, pos], [-6 6] ,'Color','r')
hold off
end

figure(111)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');   
 
subplot(2,2,1)  
hold on
plot(x_sort(:,end))
line([pos, pos], [-2 3] ,'Color','r')
hold off
title('$$x_T$$ (the last state)')
plotTickLatex2D;


subplot(2,2,2)    
hold on
plot(eps_sim_sort)
line([pos, pos], [-2 2] ,'Color','r')
hold off
title('eps\_sim ($$x_T = y^*$$ - eps\_sim)')
    plotTickLatex2D;
 
subplot(2,2,3)  
hold on
plot(eps_bar_sort)
line([pos, pos], [0 3] ,'Color','r')
hold off
title('eps\_bar (the mean of eps\_sim)')
    plotTickLatex2D;

subplot(2,2,4)  
hold on
plot(C_sim_sort)
line([pos, pos], [0 0.5] ,'Color','r')
hold off
title('C\_sim (the variance of eps\_sim)')
    plotTickLatex2D;

suptitle('Simulation Smoother: `empirical study`, NOT augmented model')

comment = {'HL vs non-HL states:','The distributions are the same,','only the variates differ'};
dim = [0.05  0.9  0.1 0.1];
% dim - four-element vector of the form [x y w h]. 
% x, y - the position; w, h - the size.
annotation('textbox',dim,'String',comment,'BackgroundColor','y');%,'FitBoxToText','on');


name = ['figures/',model,'_sim_smooth.png'];
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')