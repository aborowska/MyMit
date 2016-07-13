%% Initialisation
% clc
clear all
close all
% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s); 
addpath(genpath('include/'));

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014b)')
    v_new = 1;
else
    v_new = 0;
end

plot_on = false;
print_on = false;
save_on = true;
plot_on2 = false;

% User prompts
usr_prompt = 0;

if usr_prompt
    model = 'Select model [sv/svt]: ';
    model = input(model,'s');

    N_sim = 'Select no. of MC replications: [5/10/20] ';
    N_sim = input(N_sim);
   
    hp = 'Select no. of days-ahead for VaR estimation: [1] ';
    hp = input(hp);

    p_bar = 'Select the quantile for VaR estimation: [0.01/0.02/0.05] ';
    p_bar = input(p_bar);
else
    model = 'svt';
    N_sim = 100;
    hp = 1;
    p_bar = 0.01;
end
sim = 1;
algo = 'Direct';

%% Constants
x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

% prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -2.5*log(0.025), -log(gamma(2.5))];
% prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -0.5*log(2), -log(gamma(0.5))];
if strcmp(model, 'sv')
    prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  2.5*log(0.025), -log(gamma(2.5))];
else
    prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  2.5*log(0.025), -log(gamma(2.5)), 0.1]; % the last one is lambda for nu
end

logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
% % logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
logpdf_invgamma = @(x) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(x) - 0.025./x;
% % logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) -0.5*log(x) - 0.5*x;
if strcmp(model, 'svt')
    logpdf_exp = @(x) log(prior_const(1,5)) - prior_const(1,5)*(x - 2);  
end

%% Data
y = csvread('IBM_ret.csv');
y = 100*y;
% fprintf('IBM logreturns kurtosis: %6.4f, and skewneN_sim: %6.4f.\n', kurtosis(data), skewneN_sim(data));
% http://faculty.chicagobooth.edu/nicholas.polson/research/papers/jpr2.pdf

% QERMit ARCH
% y = csvread('GSPC_ret.csv');
% y = 100*y;
% ind_arch = find(y<=-5.5, 1, 'last' );
% y = y(1:ind_arch,1);
% y = y - mean(y);

T = size(y,1);

% SV_plot0;

%% Initialisation
EMitISEM_Control
cont.resmpl_on = false;
cont.mit.CV_tol = 0.1;
 

M = 10000;

par_NAIS_init.b = zeros(T,1);
par_NAIS_init.C = ones(T,1); 

VaR_direct = zeros(N_sim,1);
ES_direct = zeros(N_sim,1);
accept_direct = zeros(N_sim,1);
 
  
%% SML 
% load('SML_ibm.mat', 'par_SV_opt', 'heN_sim_SV_opt') 
% theta = [c, phi, sigma2, nu]
if strcmp(model,'sv')
    mu_init = [0.5, 0.98, 0.15^2];
    load('SML_ibm.mat', 'par_SV_opt', 'V_SV_corr_opt') 
%     load('SML_arch.mat', 'par_SV_opt', 'V_SV_corr_opt') 
else
    mu_init = [0.5, 0.98, 0.15^2, 7];
    load('SMLt_ibm.mat', 'par_SV_opt', 'V_SV_corr_opt') 
end

d = size(mu_init,2);

% Sigma = inv(hess_sim_SV_opt);
Sigma = V_SV_corr_opt;
Sigma = reshape(Sigma,1,d^2);
mit_init.mu = par_SV_opt;
mit_init.Sigma = Sigma;
mit_init.p = cont.mit.pnc;
mit_init.df = cont.mit.dfnc;
    
if strcmp(model,'sv')
    kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);
else
    kernel = @(a) posterior_svt(y, a, par_NAIS_init, prior_const, cont.nais);
end

mit = mit_init;
% if strcmp(model,'svt')
%     mit.mu(4) = 8;
%     mit.Sigma(16) = 9;
% end
%% Generate set of draws of theta using independence MH with naive candiate
% (based on the SML estimates) 
% simulate returns based on the draw of theta and future return paths     
% compute VaR_direct

for sim = 1:N_sim    
    fprintf('NSE sim = %i.\n', sim);
        
    [theta, x, lnw, lnk, lng_y, lnw_x, ~, accept_direct(sim,1)]= EMit_MH(M+1000, d, kernel, mit, GamMat, true);
    fprintf('MH acceptance rate: %6.4f (%s, %s). \n', accept_direct(sim,1), model, algo);

    theta = theta(1001:M+1000,:);
    x = x(1001:M+1000,:);
    lnw = lnw(1001:M+1000,:);
    lnk = lnk(1001:M+1000,:);

    if (sim == N_sim)
        SV_autocorr;
    end
    % High loss, 1 days horizon
    % approximate the high loss distribution of (theta,eps*,eta*) where 
    % eps*={eps_T+1,...,eps_T+hp}
    % eta*={eta_T+1,...,eta_T+hp}

    c = theta(:,1);
    phi = theta(:,2);
    sigma2 = theta(:,3);
    
    eta_h1 = randn(M,1);
    
    if strcmp(model,'sv')
        eps_h1 = randn(M,1);
        x_h1 = c + phi.*(x(:,end) - c) + sqrt(sigma2).*eta_h1;
        y_h1 = exp(0.5*x_h1).*eps_h1;
    else
        nu = theta(:,4);
        rho = (nu-2)./nu;
        eps_h1 = trnd(repmat(nu,1,hp));

        x_h1 = c + phi.*(x(:,end) - c) + sqrt(sigma2).*eta_h1;
        y_h1 = sqrt(rho).*exp(0.5*x_h1).*eps_h1;
    end
    
    ind_real = (imag(y_h1)==0);
    M_real = sum(ind_real); 
    y_h1 = y_h1(ind_real,:);
        
    [PL_h1, ind] = sort(fn_PL(y_h1));
   
    VaR_direct(sim,1) = PL_h1(round(p_bar*M_real));
    ES_direct(sim,1) = mean(PL_h1(round(1:p_bar*M)));   
    fprintf('p_bar = %4.2f, VaR_direct  = %4.5f. \n', p_bar, VaR_direct(sim,1))
    fprintf('p_bar = %4.2f, VaR_direct  = %4.5f. \n', p_bar, mean(VaR_direct(VaR_direct<0,1)))
    fprintf('p_bar = %4.2f, NSE VaR_direct  = %4.5f. \n', p_bar, std(VaR_direct(VaR_direct<0,1)))
end

 
if strcmp(model,'sv')
    save(['results/sv_VaR_direct_',int2str(M),'.mat'],'mit','accept_direct','theta', 'x', 'lnw', 'lnk','p_bar','M','N_sim','VaR_direct','ES_direct','ind','ind_real');
else
    save(['results/svt_VaR_direct_',int2str(M),'.mat'],'mit','accept_direct','theta', 'x', 'lnw', 'lnk','p_bar','M','N_sim','VaR_direct','ES_direct','ind','ind_real');
end


%%
mean_VaR_direct = mean(VaR_direct);
mean_ES_direct = mean(ES_direct);

NSE_VaR_direct = std(VaR_direct);
NSE_ES_direct = std(ES_direct);
NSE_ES_IS = std(ES_direct);

fprintf('(%s) 100*%4.2f%% VaR direct (mean) estimate: %6.4f. \n', model, p_bar, mean_VaR_direct);
fprintf('(%s) NSE VaR direct: %6.4f. \n', model, NSE_VaR_direct);
fprintf('(%s) VaR direct: [%6.4f, %6.4f]. \n \n', model, mean_VaR_direct - NSE_VaR_direct, mean_VaR_direct + NSE_VaR_direct);
 

fprintf('(%s) 100*%4.2f%% ES direct (mean) estimate: %6.4f. \n', model, p_bar, mean_ES_direct);
fprintf('(%s) NSE ES direct: %6.4f. \n', model, NSE_ES_direct);
fprintf('(%s) ES direct: [%6.4f, %6.4f]. \n \n', model, mean_ES_direct - NSE_ES_direct, mean_ES_direct + NSE_ES_direct);

 

if plot_on2
    figure(590+100*p_bar)
%         set(gcf, 'visible', 'off');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);   
    set(gcf,'defaulttextinterpreter','latex');
    boxplot(VaR_direct,'labels',{'VaR_direct'})        
    title(['100*', num2str(p_bar),'\% VaR estimates: prelim and IS (',strrep(model,'_','\_'),', ',algo,', M = ',num2str(M),', N\_sim = ', num2str(N_sim),').'])  
    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    if print_on
        name = ['figures/',model,'_',algo,'_', num2str(p_bar),'_VaR_box_',num2str(M),'.png'];
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%

    figure(5900+100*p_bar)
%         set(gcf, 'visible', 'off');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);   
    set(gcf,'defaulttextinterpreter','latex');
    hold on; 
    bar(VaR_IS,'FaceColor',[0 0.4470 0.7410], 'EdgeColor','w'); 
    plot(0:(N_sim+1), (mean_VaR_prelim - NSE_VaR_prelim)*ones(N_sim+2,1),'r--'); 
    plot(0:(N_sim+1), (mean_VaR_prelim + NSE_VaR_prelim)*ones(N_sim+2,1),'r--'); 
    plot(0:(N_sim+1), mean_VaR_prelim*ones(N_sim+2,1),'r'); 
    hold off;
    title(['(',model,',',algo,' M = ',num2str(M),') ','100*', num2str(p_bar),'\% VaR IS estimates and the mean VaR prelim (+/- NSE VaR prelim).'])

    if v_new
        set(gca,'TickLabelInterpreter','latex')
    else
        plotTickLatex2D;
    end
    if print_on
        name = ['figures/',model,'_',algo,'_', num2str(p_bar),'_VaR_bar_',num2str(M),'.png'];
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    figure(777)
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
    set(gcf,'defaulttextinterpreter','latex');
    hold on
    plot(PL_h1,'b')
    pos =  max(find(PL_h1<=mean(VaR_direct)));
    scatter(pos, mean(VaR_direct),'MarkerEdgeColor','green','MarkerFaceColor','green')
    hold off
    plotTickLatex2D;

    title(['Sorted future profit/losses values $$PL(y_{T+1}^{(i)})$$. Model: $$',strrep(model,'_','\_'),'$$.'])
    if print_on
        if strcmp(model,'sv_x')
            name = 'figures/sv_x_predict_direct.png';
        elseif strcmp(model,'sv')
            name = 'figures/sv_predict_direct.png';
        else
            name = 'figures/svt_predict_direct.png';
        end
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(name,'-dpng','-r0')
    end    

end

if save_on
    gen_out2;
end