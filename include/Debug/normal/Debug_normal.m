clear all
close all
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

p_bar = 0.01;
M = 10000;

MitISEM_Control
cont.mit.dfnc = 5;
cont.mit.N = 10000;
cont2 = cont;
N_sim = 10;

kernel_init =  @(aa) 0.5*(log(2*pi) + aa.^2);
kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
mu_init = 1;

% options = optimset('Display','iter');
% [x0,val,~,~,~,H] = fminunc(kernel_init,mu_init,options);
% fminsearch(kernel_init,mu_init)

mit_init.mu = 0;
mit_init.Sigma = 1;
mit_init.df = 1;
mit_init.p = 1;
[mit1, summary1] = MitISEM_new(mit_init, kernel, mu_init, cont, GamMat);

figure(1)
xx = -3:0.01:3;
Mit1 = dmvgt(xx', mit1, false, GamMat);    
hold on
plot(xx,exp(kernel(xx)), 'g')
hold off

[draw, lnk1, ~] = fn_rmvgt_robust(M, mit1, kernel);
[draw1, accept] = Mit_MH(M+1000, kernel, mit1, GamMat);
draw1 = draw1(1001:M+1000);
[draw1_sort, ind] = sort(draw1);
VaR_prelim = draw1_sort(p_bar*M); 
ES_prelim = mean(draw1_sort(1:p_bar*M));   
VaR_true = norminv(p_bar);
ES_true = -normpdf(VaR_true)/p_bar;
  
VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);
% for sim =1:N_sim
kernel_init = @(aa) - tail_normal(aa,VaR_prelim);
kernel = @(aa) tail_normal(aa,VaR_prelim);
mu_hl = draw1_sort(max(find(draw1_sort < VaR_prelim)),:);

[mit2, summary2] = MitISEM_new(kernel_init, kernel, mu_hl, cont2, GamMat);

figure(2)
xx = -5:0.01:5;
Mit2 =  dmvgt(xx', mit2, false, GamMat);  
plot(xx,Mit2)
hold on
plot(xx,exp(kernel(xx'))/p_bar, 'g')
hold off

for sim =1:N_sim
    kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
    [draw1, lnk1, ~] = fn_rmvgt_robust(M, mit1, kernel); 
    kernel = @(aa) tail_normal(aa,VaR_prelim);
    [draw2, lnk2, ~] = fn_rmvgt_robust(M, mit2, kernel);
    draw_opt = [draw1; draw2];
    lnk_opt = [lnk1; lnk2];
    
    exp_lnd1 = 0.5*dmvgt(draw_opt,mit1,false, GamMat);
    exp_lnd2 = 0.5*dmvgt(draw_opt,mit2,false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd); % take log to comute the importance weights in fn_ISwgts
    w_opt =  fn_ISwgts(lnk_opt, lnd_opt, false); % false - not normalised --> will be in fn_PL function

    hl_w  = sum( w_opt( draw_opt(:,1)<VaR_prelim,:)/sum(w_opt) );    

    dens = struct('y',draw_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1); %  computed for the case when (L == 1) 
    VaR_IS(sim) = IS_estim(1,1);
    ES_IS(sim) = IS_estim(1,2);
end  

plot(VaR_IS)
hold on
plot(ones(N_sim,1)*VaR_prelim,'r')
plot(ones(N_sim,1)*VaR_true,'g')
hold off

