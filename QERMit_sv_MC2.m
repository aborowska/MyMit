N_sim = 10;
M = 2000;
VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONTE CARLO NSE AND RNE ESTIMATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sim = 1:N_sim    
    fprintf('NSE sim = %i.\n', sim);

    theta1 = rmvgt2(M, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
    eta_h1_1 = randn(M,1);
    if strcmp(model,'sv')
        eps_h1_1 = randn(M,1);
    else
        nu1 = theta1(:,4);
        rho1 = (nu1-2)./nu1;
        eps_h1_1 = trnd(repmat(nu1,1,hp));
    end      
    theta1 = [theta1, eta_h1_1, eps_h1_1];
    theta2 = rmvgt2(M, mit2.mu, mit2.Sigma, mit2.df, mit2.p); 
    theta_opt = [theta1; theta2];
    lnk_opt = zeros(2*M,1);
    x_end = zeros(2*M,1);
    for ii = 1:(2*M)
        [par_NAIS, x_smooth]= NAIS_param(par_NAIS_init, y, theta_opt(ii,1:3), cont.nais);
        [x, lng_y, lnw] = NAIS_loglik(y, theta_opt(ii,1:3), par_NAIS, cont.nais);
        x_end(ii,1) = x(end);
        prior = prior_sv(theta_opt(ii,1:3),prior_const);
        prior_eta =  prior_const(1,1) - 0.5*(theta_opt(ii,4)).^2;
        prior_eps =  prior_const(1,1) - 0.5*(theta_opt(ii,5)).^2;
        lnk_opt(ii,1) = lng_y + lnw + prior + prior_eta + prior_eps;
    end
    exp_lnd1 = 0.5*exp(2*prior_const(1,1) - 0.5*(theta_opt(:,4)).^2 - 0.5*(theta_opt(:,5)).^2 + dmvgt(theta_opt(:,1:3), mit1, true, GamMat));
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
        x_opt_h1 = c_opt + phi_opt.*(x_end - c_opt) + sqrt(sigma2_opt).*eta_opt;
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
    PL_opt(imag(PL_opt) ~= 0) = -Inf;
    [PL_opt_h1, ind] = sort(PL_opt);
    w_opt_h1 = w_opt(ind)/sum(w_opt);
    cum_w = cumsum(w_opt_h1);
    ind_var = sum(cum_w<=p_bar);
    VaR_IS(sim,1) = PL_opt_h1(ind_var);
end

boxplot([VaR_prelim_MC, VaR_IS],'labels',{'VaR_prelim MC','VaR_IS'})        
   