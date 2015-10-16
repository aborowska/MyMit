N_sim = 10;
cont.mit.N = 2000;
N = cont.mit.N;
    
VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONTE CARLO NSE AND RNE ESTIMATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sim = 1:N_sim
    fprintf('\n')
    fprintf('NSE sim = %i.\n', sim);
    fprintf('\n')
    
    %%% MIT1 %%%
    if strcmp(model,'sv')
        kernel_prior = @(a) prior_sv(a, prior_const); 
        kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);
    else
        kernel_prior = @(a) prior_svt(a, prior_const); 
        kernel = @(a) posterior_svt(y, a, par_NAIS_init, prior_const, cont.nais);
    end
    [theta1, x1, lnw1, lnk1, lng_y1, lnw_x1, x_smooth1, ~] = EMit_MH(y, N, kernel_prior, kernel, mit1, GamMat, false);

    eta_h1_1 = randn(N,1);

    if strcmp(model,'sv')
        eps_h1_1 = randn(N,1);
    else
        nu1 = theta1(:,4);
        rho1 = (nu1-2)./nu1;
        eps_h1_1 = trnd(repmat(nu1,1,hp));
    end    
    
    %%% MIT2 %%% 
    if strcmp(model,'sv')
        kernel_prior = @(a) prior_sv_hl(a, prior_const); 
        kernel = @(a) posterior_sv_hl(y, a, VaR_prelim, par_NAIS_init, prior_const, cont.nais); 
    else
        kernel_prior = @(a) prior_svt_hl(a, prior_const); 
        kernel = @(a) posterior_svt_hl(y, a, VaR_prelim, par_NAIS_init, prior_const, cont.nais); 
    end

    [theta2, x2, lnw2, lnk2, lng_y2, lnw_x2, x_smooth2, ~] = EMit_MH(y, N, kernel_prior, kernel, mit2, GamMat, false);

    %%% OPT %%%
    theta_opt = [theta1, eta_h1_1, eps_h1_1; theta2];
    x_opt_end = [x1(:,end); x2(:,end)];

    if strcmp(model,'sv')
        lnk_opt = lnk1 + 2*prior_const(1,1) - 0.5*(eta_h1_1).^2 - 0.5*(eps_h1_1).^2;
    else
        lnk_opt = lnk1 + prior_const(1,1) - 0.5*(eta_h1_1).^2 + duvt(eps_h1_1, theta1(:,4), hp, true);
    end
    lnk_opt = [lnk_opt; lnk2];

% lng_y_opt = [lng_y1; lng_y2];
% lnw_x_opt = [lnw_x1; lnw_x2];
    %% IS weights
    if strcmp(model,'sv')
        exp_lnd1 = 0.5*normpdf(theta_opt(:,4)).*normpdf(theta_opt(:,5)).*dmvgt(theta_opt(:,1:3), mit1, false, GamMat);
    else
        exp_lnd1 = 0.5*normpdf(theta_opt(:,5)).*duvt(theta_opt(:,6), theta_opt(:,4), 1, false).*dmvgt(theta_opt(:,1:4), mit1, false, GamMat);
    end
    exp_lnd2 = 0.5*dmvgt(theta_opt, mit2, false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);

    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    c_opt = theta_opt(:,1);
    phi_opt = theta_opt(:,2);
    sigma2_opt = theta_opt(:,3);
 
    if strcmp(model,'sv')
        eta_opt = theta_opt(:,4);
        eps_opt = theta_opt(:,5);  
        x_opt_h1 = c_opt + phi_opt.*(x_opt_end - c_opt) + sqrt(sigma2_opt).*eta_opt;
        y_opt_h1 = exp(0.5*x_opt_h1).*eps_opt;
    else
        nu_opt = theta_opt(:,4);
        rho_opt = (nu_opt-2)./nu_opt;
        eta_opt = theta_opt(:,5);
        eps_opt =  theta_opt(:,6);
        x_opt_h1 = c_opt + phi_opt.*(x_opt_end - c_opt) + sqrt(sigma2_opt).*eta_opt;
        y_opt_h1 = sqrt(rho_opt).*exp(0.5*x_opt_h1).*eps_opt;
    end

    PL_opt = fn_PL(y_opt_h1);
    [PL_opt_h1, ind] = sort(PL_opt);
    w_opt_h1 = w_opt(ind)/sum(w_opt);
    cum_w = cumsum(w_opt_h1);
%  
% lnk_opt_h1 = lnk_opt(ind);
% lnd_opt_h1 = lnd_opt(ind);
% lng_y_opt_h1 = lng_y_opt(ind);
% lnw_x_opt_h1 = lnw_x_opt(ind); 
% lnprior_opt = lnk_opt - lng_y_opt - lnw_x_opt;
% lnprior_opt_h1 = lnk_opt_h1 - lng_y_opt_h1 - lnw_x_opt_h1;
% theta_opt_h1 = theta_opt(ind,:);
% x_opt = [x1; x2];
% x_opt_h1 = x_opt(ind,:);
% 
% x_smooth = [x_smooth1, x_smooth2];
% x_smooth_h1 = x_smooth(:,ind);
% 
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'defaulttextinterpreter','latex');
subplot(2,2,1); plot(cum_w); hold on; plot(0.01*ones(2*N),'r'); hold off; title('cum w');
subplot(2,2,2); plot(lnk_opt_h1);   title('lnk opt h1');
subplot(2,2,3); plot(w_opt_h1);     title('w opt h1');
subplot(2,2,4); plot(lnd_opt_h1);   title('lnd opt h1');
% 
% figure(2)
% subplot(2,2,1); plot(lnk2); title('lnk2');
% subplot(2,2,2); plot(lng_y2); title('lngy2');
% subplot(2,2,3); plot(lnd_opt(2001:end)); title('lnd2');
% subplot(2,2,4); plot(lnw_x2); title('lnwx2');
% 
% 
% figure(3)
% subplot(2,2,1); plot(lnprior_opt); title('lnprior opt');
% subplot(2,2,2); plot(lng_y_opt); title('lngy opt');
% subplot(2,2,3); plot(lnd_opt); title('lnd opt');
% subplot(2,2,4); plot(lnw_x_opt); title('lnwx opt');
% 
% 
% [~,ind10]=sort(lnk_opt_h1);
% ind10 = ind10(1:10);
% theta10=theta_opt_h1(ind10,:);

% subplot(2,3,1)
% plot(w_opt_h1)
% title('w_opt_h1')
% subplot(2,3,2)
% plot(lnk_opt_h1)
% title('lnk_opt_h1')
% subplot(2,3,3)
% plot(lnd_opt_h1)
% title('lnd_opt_h1')
% subplot(2,3,4)
% plot(lng_y_opt_h1)
% title('lng_y_opt_h1')
% subplot(2,3,5)
% plot(lnw_x_opt_h1)
% title('lnw_x_opt_h1')
% subplot(2,3,6)
% plot(lnprior_opt_h1)
% title('lnprior_opt_h1')
    

    dens = struct('y',y_opt_h1,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);

    SV_plot3;

end

mean_VaR_IS = mean(VaR_IS(VaR_IS<0));
mean_ES_IS = mean(ES_IS(ES_IS<0));

NSE_VaR_IS = std(VaR_IS(VaR_IS<0));              
NSE_ES_IS = std(ES_IS(ES_IS<0));                
           

fprintf('IS VAR (mean) estimate: %6.4f. \n',mean_VaR_IS);
fprintf('IS ES (mean) estimate: %6.4f. \n',mean_ES_IS);

fprintf('IS NSE VaR estimate: %6.4f. \n',NSE_VaR_IS);
fprintf('IS NSE ES estimate: %6.4f. \n',NSE_ES_IS);

clear GamMat x_gam fig h f xi xx

if strcmp(model,'sv')
    save(['results/sv_nse.mat']);
else
    save(['results/svt_nse.mat']);
end