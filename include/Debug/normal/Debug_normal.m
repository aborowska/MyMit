%%
clear all
close all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

p_bar = 0.05;
M = 10000;

MitISEM_Control
cont.mit.dfnc = 5;
cont.mit.N = 10000;
cont2 = cont;
N_sim = 1000;

kernel_init =  @(aa) 0.5*(log(2*pi) + aa.^2);
kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
mu_init = 1;

VaR_true = norminv(p_bar);
ES_true = -normpdf(VaR_true)/p_bar;

%% mit1
mit_init.mu = 0;
mit_init.Sigma = 1;
mit_init.df = 1;
mit_init.p = 1;

[mit1, summary1] = MitISEM_new(mit_init, kernel, mu_init, cont, GamMat);

%% figure(1)
figure(1)
xx = -4:0.01:4;
Mit1 = dmvgt(xx', mit1, false, GamMat);  
plot(xx,Mit1)
hold on
plot(xx,exp(kernel(xx)), 'g')
hold off
legend('mit1','kernel')
% name = 'include/debug/normal/normal_posterior.png';
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% print(name,'-dpng','-r0')

%% VaR_prelim
VaR_prelim_MC = zeros(N_sim,1);
ES_prelim_MC = zeros(N_sim,1);

resampl_on = false;
for sim = 1:N_sim
    [draw, ~, ~] = fn_rmvgt_robust(M, mit1, kernel, resampl_on);
    [draw1, accept] = Mit_MH(M+1000, kernel, mit1, GamMat);
    draw1 = draw1(1001:M+1000);
    [draw1_sort, ind] = sort(draw1);
    VaR_prelim_MC(sim,1) = draw1_sort(p_bar*M); 
    ES_prelim_MC(sim) = mean(draw1_sort(1:p_bar*M));   
end

%% mit2
VaR_prelim = mean(VaR_prelim_MC);
ES_prelim = mean(ES_prelim_MC);

% for sim =1:N_sim
kernel_init = @(aa) - tail_normal(aa,VaR_prelim);
kernel = @(aa) tail_normal(aa,VaR_prelim);
mu_hl = draw1_sort(max(find(draw1_sort < VaR_prelim)),:);

[mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);

%% figure(2)
figure(2)
xx = -5:0.01:5;
Mit2 =  dmvgt(xx', mit2, false, GamMat);  
plot(xx,Mit2)
hold on
plot(xx,exp(kernel(xx'))/p_bar, 'g')
hold off
legend('mit2','kernel')
name = 'include/debug/normal/normal_tail.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

%% VaR_IS
VaR_IS = zeros(N_sim,1);
VaR_IS_kernel = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);

M1 = zeros(N_sim,1);
M2 = zeros(N_sim,1);
M3 = zeros(N_sim,1);
M4 = zeros(N_sim,1);

M1_draw1 = zeros(N_sim,1);
M2_draw1 = zeros(N_sim,1);
M3_draw1 = zeros(N_sim,1);
M4_draw1 = zeros(N_sim,1);

M1_draw2 = zeros(N_sim,1);
M2_draw2 = zeros(N_sim,1);
M3_draw2 = zeros(N_sim,1);
M4_draw2 = zeros(N_sim,1);


resampl_on = false;
for sim = 1:N_sim
    kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
    [draw1, lnk1, ~] = fn_rmvgt_robust(M, mit1, kernel, resampl_on); 
    M1_draw1(sim,1) = mean(draw1); % sum(draw1.^1)/length(draw1);
    M2_draw1(sim,1) = var(draw1); % sum((draw1.^2-mean(draw1)))/length(draw1);
    M3_draw1(sim,1) = skewness(draw1); % sum(draw1.^3)/length(draw1);
    M4_draw1(sim,1) = kurtosis(draw1); % sum(draw1.^4)/length(draw1);
    
%     kernel = @(aa) tail_normal(aa,VaR_true);
    [draw2, lnk2, ~] = fn_rmvgt_robust(M, mit2, kernel, resampl_on);
    
%     kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
%     lnk22 = kernel(draw2);
    
    M1_draw2(sim,1) = mean(draw2); % sum(draw2.^1)/length(draw1);
    M2_draw2(sim,1) = var(draw2); % sum(draw2.^2)/length(draw1);
    M3_draw2(sim,1) = skewness(draw2); % sum(draw2.^3)/length(draw1);
    M4_draw2(sim,1) = kurtosis(draw2); % sum(draw2.^4)/length(draw1);
    
    draw_opt = [draw1; draw2];
    [draw_sort, ind] = sort(draw_opt); 

    lnk_opt = [lnk1; lnk2];
%     draw_opt = 0.5*draw1 + 0.5*draw2;
%     kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
%     lnk_opt = kernel(draw_opt);

    lnd_opt = log(0.5*dmvgt(draw_opt,mit1,false, GamMat) + 0.5*dmvgt(draw_opt,mit2,false, GamMat));

%     exp_lnd1 = 0.5*dmvgt(draw_opt,mit1,false, GamMat);
%     exp_lnd2 = 0.5*dmvgt(draw_opt,mit2,false, GamMat);
%     exp_lnd = exp_lnd1 + exp_lnd2;
%     lnd_opt = log(exp_lnd); % take log to comute the importance weights in fn_ISwgts
      
%     w_opt =  fn_ISwgts(lnk_opt, lnd_opt, false); % false - not normalised --> will be in fn_PL function
w_opt = lnk_opt-lnd_opt;
w_opt = exp(w_opt);    
w = w_opt(ind);
cum_w = cumsum(w);
M1(sim,1) = sum((draw_opt).*w_opt)/length(draw_opt);
M2(sim,1) = sum((draw_opt.*w_opt - M1(sim,1)).^2)/length(draw_opt);
M3(sim,1) = sum((draw_opt.^3).*w_opt)/length(draw_opt);
M4(sim,1) = sum((draw_opt.^4).*w_opt)/length(draw_opt);

%     w = w_opt(ind,:);
%     cum_w = cumsum(w);
ind_kernel = sum(cumsum(w/sum(w)) < p_bar);
    ind_var = min(find(cum_w > length(draw_opt)*p_bar))-1; 
ind_exact = sum(cumsum(w/length(draw_opt)) < p_bar);    
% fprintf('ind_var: %i, ind_kernel: %i, ind_exact: %i\n', ind_var, ind_kernel, ind_exact);
%     if isempty(ind_var)
%         ind_var = length(y)-1;
%     end
%     if (ind_var == 0)
%         ind_var = 1;
%     end
    VaR_IS(sim,1) = draw_sort(ind_var); % intrapolate
    VaR_IS_kernel(sim,1) = draw_sort(ind_kernel); % intrapolate

    %     VaR_IS(sim,1) = (draw_sort(ind_var+1) + draw_sort(ind_var))/2; % intrapolate
    ES_IS(sim,1) = sum((w(1:ind_var)/sum(w(1:ind_var))).*draw_sort(1:ind_var));
end

%% figure(3)
figure(3)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
plot(VaR_IS)
hold on
plot(VaR_IS_kernel,'c')
plot(VaR_prelim_MC,'r')
plot(ones(N_sim,1)*VaR_true,'g')
plot(ones(N_sim,1)*VaR_prelim,'k')
hold off
legend('VaR\_IS','VaR\_IS\_kernel','VaR\_prelim\_MC','VaR\_true','VaR\_used')
name = 'include/debug/normal/normal_VaR.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')


figure(4)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
plot(VaR_IS)
hold on
plot(VaR_IS_kernel,'c')
% plot(VaR_prelim_MC,'r')
plot(ones(N_sim,1)*VaR_true,'g')
plot(ones(N_sim,1)*VaR_prelim,'k')
hold off
legend('VaR\_IS','VaR\_IS\_kernel', 'VaR\_true','VaR\_used')
name = 'include/debug/normal/normal_VaR2.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')

figure(5)
boxplot([VaR_prelim_MC, VaR_IS],'labels',{'VaR_prelim MC','VaR_IS'})        
name = 'include/debug/normal/normal_VaR_box.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')