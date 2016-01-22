
f_pl = @(aa) 100*(exp(aa/100) - 1); 

ind_y = (imag(y_opt_h1)==0); 
% y_save = y_opt_h1;
y_opt_h1 = y_opt_h1(ind_y);
% w_save = w_opt;
% w_opt = w_opt(ind_y);
PL_opt_h1 = f_pl(sum(y_opt_h1,2));
[PL_opt_h1, ind] = sort(PL_opt_h1); 
% w_opt = w_opt(ind);
y_opt_h1 = y_opt_h1(ind);
% plot(y_opt_h1,w_opt)

figure(777)
set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
set(gcf,'defaulttextinterpreter','latex');
hold on
plot(PL_opt_h1,'b')
pos =  max(find(PL_opt_h1<= mean(VaR_IS(VaR_IS<0))));
scatter(pos,  mean(VaR_IS(VaR_IS<0)),'MarkerEdgeColor','red','MarkerFaceColor','red')
pos =  max(find(PL_opt_h1<=VaR_prelim));
scatter(pos, VaR_prelim,'MarkerEdgeColor','green','MarkerFaceColor','green')
hold off
plotTickLatex2D;

title(['Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$. Model: $$',strrep(model,'_','\_'),'$$.'])
if strcmp(model,'sv_x')
    name = 'figures/sv_x_predict.png';
elseif strcmp(model,'sv')
    name = 'figures/sv_predict.png';
else
	name = 'figures/sv_t_predict.png';
end
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

            
ind_x = ~isnan(x);
x_opt_h1 = x_opt_h1(ind_y,:);
x = x(ind_y,:);
eps_bar = eps_bar(ind_y,:);
eps_sim = eps_sim(ind_y,:);
C_sim = C_sim(ind_y,:);

x_T = x(ind,:);
eps_bar_T = eps_bar(ind,1);
eps_sim_T = eps_sim(ind,1);
C_sim_T = C_sim(ind,1);

pos =  max(find(PL_opt_h1<=VaR_IS(sim,1)));

x_T = x_T(1:pos,:);
eps_bar_T = eps_bar_T(1:pos,:);
eps_sim_T = eps_sim_T(1:pos,:);
C_sim_T = C_sim_T(1:pos,:);


figure(10010)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');
subplot(2,5,1)
    hist(x_T)
    title('x\_T all')
subplot(2,5,6)
    hist(x_T_hl_init)
    title('x\_T high loss')
subplot(2,5,2)
    hist(eps_bar)
    title('eps\_bar all')
subplot(2,5,7)
    hist(eps_bar_hl_init)
    title('eps\_bar high loss')
subplot(2,5,3)
    hist(eps_sim)
    title('eps\_sim all')
subplot(2,5,8)
    hist(eps_sim_hl_init)
    title('eps\_bar high loss')
subplot(2,5,4)
    hist(C_T)
    title('C all')
subplot(2,5,9)
    hist(C_T_hl_init)
    title('C high loss')
subplot(2,5,5)
    hist(RND)
    title('RND all')
subplot(2,5,10)
    hist(RND_hl_init)
    title('RND high loss')

fprintf('x_T all >> mean: %6.4f, std: %6.4f.\n',mean(x_T(~isnan(x_T))), std(x_T(~isnan(x_T))))    
fprintf('x_T h_l >> mean: %6.4f, std: %6.4f.\n\n',mean(x_T_hl_init(~isnan(x_T_hl_init))), std(x_T_hl_init(~isnan(x_T_hl_init))))    

fprintf('eps_bar all >> mean: %6.4f, std: %6.4f.\n',mean(eps_bar(~isnan(eps_bar))), std(eps_bar(~isnan(eps_bar))))    
fprintf('eps_bar h_l >> mean: %6.4f, std: %6.4f.\n\n',mean(eps_bar_hl_init(~isnan(eps_bar_hl_init))), std(eps_bar_hl_init(~isnan(eps_bar_hl_init))))    

fprintf('eps_sim all >> mean: %6.4f, std: %6.4f.\n',mean(eps_sim(~isnan(eps_sim))), std(eps_sim(~isnan(eps_sim))))    
fprintf('eps_sim h_l >> mean: %6.4f, std: %6.4f.\n\n',mean(eps_sim_hl_init(~isnan(eps_sim_hl_init))), std(eps_sim_hl_init(~isnan(eps_sim_hl_init))))    

fprintf('C_T all >> mean: %6.4f, std: %6.4f.\n',mean(C_T(~isnan(C_T))), std(C_T(~isnan(C_T))))    
fprintf('C_T h_l >> mean: %6.4f, std: %6.4f.\n\n',mean(C_T_hl_init(~isnan(C_T_hl_init))), std(C_T_hl_init(~isnan(C_T_hl_init))))    

fprintf('RND all >> mean: %6.4f, std: %6.4f.\n',mean(RND(~isnan(RND))), std(RND(~isnan(RND))))    
fprintf('RND h_l >> mean: %6.4f, std: %6.4f.\n\n',mean(RND_hl_init(~isnan(RND_hl_init))), std(RND_hl_init(~isnan(RND_hl_init))))    

name = 'figures/sv_x_hist.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')
