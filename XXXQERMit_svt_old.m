 %% Initialisation
clear all
addpath(genpath('include/'));

addpath('results/');

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);
prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -2.5*log(0.025), -log(gamma(2.5))];

model = 'svt';
y = csvread('IBM_ret.csv');
y = 100*y;

% fprintf('IBM logreturns kurtosis: %6.4f, and skewness: %6.4f.\n', kurtosis(data), skewness(data));

% http://faculty.chicagobooth.edu/nicholas.polson/research/papers/jpr2.pdf
T = size(y,1);

plot_on = false;
print_on = false;

if plot_on
	figure(1)
	set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);
	set(gcf,'defaulttextinterpreter','latex');
    xx = linspace(2007,2012,T);
    plot(xx, y)
    set(gca,'XTickLabel',{'2007',' ','2008',' ','2009',' ','2010',' ','2011',' ','2012'})
	title('IBM log-returns $$y$$')
%    plotTickLatex2D;
	set(gca,'TickLabelInterpreter','latex')
	if print_on
		name = 'figures/sv_data.png';
		fig = gcf;
		fig.PaperPositionMode = 'auto';
		print(name,'-dpng','-r0')
	end
end

EMitISEM_Control
cont.mit.dfnc = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1a.:
% theta = [c, phi, sigma2, nu]
% mu_init = [0.5, 0.98, 0.15^2, 7];
mu_init = [0.5, 0.98, 0.15^2, 7];

d = size(mu_init,2);
% SML ... <<<<<<<<<<<<<<<<<<<<<<
load('SMLt_ibm.mat', 'par_SV_opt', 'hess_SV_opt') % old hessina = hessina correcte by the jacobian

Sigma = inv(hess_SV_opt);
Sigma(1,1) = 1;
% Sigma(1:3,1:3) = Sigma(1:3,1:3)/10;
Sigma(4,:) = Sigma(4,:)/20;
Sigma(:,4) = Sigma(:,4)/20;
Sigma = reshape(Sigma,1,d^2);

mit_init.mu = par_SV_opt;
mit_init.Sigma = Sigma;
mit_init.p = 1;
mit_init.df = 5;
mit = mit_init;
kernel = @(a) prior_svt(a,true,prior_const);
[mit1, theta1, x1, b1, C1, summary1]  = EMitISEM(mit_init, kernel, mu_init, cont, GamMat, prior_const);
clear GamMat
save('results/sv_mitisem.mat');
x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);

load('sv_mitisem_dom_all.mat', 'mit_adapt');

%% QERMit 1b.:
% generate set opf draws of theta using independence MH with
% candiate from MitISEM; then simulate returns based on the draw of theta 
[theta, accept] = Mit_MH(M+1000, kernel, mit1, GamMat);
fprintf('(MitISEM) MH acceptance rate: %6.4f. \n',accept);
save('results/sv_mitisem_theta.mat','theta', 'accept');

theta = theta(1001:M+1000,:);

%% High loss, 10 days horizon
% approximate the high loss distribution of (theta,eps*) where eps*={eps_T+1,...,eps_T+hp}
p_bar = 0.01;
hp = 1; % prediction horizon 

h_T = volatility_t_garch(theta, data, S);
[y_hp, eps_hp] = predict_t_garch(theta, y_T, S, h_T, hp);
% get the preliminary 10-day-ahead 99% VaR estimate as the 100th of the ascendingly sorted percentage loss values
[PL_hp, ind] = sort(fn_PL(y_hp));
eps_hl = eps_hp(ind,:); 
eps_hl = eps_hl(1:p_bar*M,:);
theta_hl = theta(ind,:);
theta_hl = theta_hl(1:p_bar*M,:);

VaR_prelim = PL_hp(p_bar*M);
fprintf('hp = %i, y_T = %4.2f, VaR_prelim = %4.5f. \n', hp, y_T, VaR_prelim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 1c.:
% get mit approximation of the conditional joint density of
% parameters and future returns given the returns are below VaR_prelim
% approximation of the joint high loss distribution
% here: not future returns but future disturbances  (varepsilons)

%% High loss density approximation
L = true;
kernel_init = @(a) - posterior_t_garch_hl(a, data, S, VaR_prelim, L, hyper, GamMat);
kernel = @(a) posterior_t_garch_hl(a, data, S, VaR_prelim, L, hyper, GamMat);
% theta = [alpha, beta, mu, nu, eps_T1, ..., eps_Thp] <-- size 14 
mu_init_hl = [0.07, 0.93, 0.05, 8.5, -1.5*ones(1, hp)];
% mu_init_hl = [0.0628    0.9311    0.0491    8.2365   -2.9356];
% mu_hl = median([theta_hl, eps_hl],1);
% kernel(mu_hl)

%% IS estimation of the initial mixture component 
draw_hl = [theta_hl, eps_hl];
% lnk = fn_lkereval(kernel, dupa);
lnk = kernel(draw_hl);
ind = find(lnk~=-Inf);
lnk = lnk(ind,:);
draw_hl = draw_hl(ind,:);
% lnd = dmvgt(theta,mit_init,true);
% w = fn_ISwgts(lnk, lnd, norm);
w = lnk/sum(lnk);
[mu_hl, Sigma_hl] = fn_muSigma(draw_hl, w);

mit_hl.mu = mu_hl;
mit_hl.Sigma = Sigma_hl;
mit_hl.df = 5;
mit_hl.p = 1;

cont2 = cont;  
cont2.mit.Hmax = 2;
[mit2, summary2] = MitISEM(mit_hl, kernel, mu_init_hl, cont2, GamMat);
%[mit2, summary] = MitISEM(kernel_init, kernel, mu_hl, cont, GamMat);
save(['results/sv_mitisem_hl_',num2str(hp),'.mat'],'mit2','summary2');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QERMit 2:
% use the mixture 0.5*mit1 + 0.5*mit2 as the importance density
% to estiamte VaR and ES for theta and y (or alpha in eps)

plot_on = false;
print_on = false;
VaR_IS = zeros(N_sim,1);
ES_IS = zeros(N_sim,1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MONTE CARLO NSE AND RNE ESTIMATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sim = 1:N_sim
    fprintf('NSE sim = %i.\n', sim);
    kernel = @(a) posterior_t_garch(a, data, S, L, hyper, GamMat);
    [draw1, lnk1, ind_red1] = fn_rmvgt_robust(M/2, mit1, kernel);
    eps1 = zeros(M/2, hp);
    for hh = 1:hp
       eps1(:,hh) = trnd(draw1(:,4)); % ERRORS ARE iid T!!
    end
    draw1_eps1 = [draw1(1:M/2,:), eps1];
    
    kernel = @(a) posterior_t_garch_hl(a, data, S, VaR_prelim, L, hyper, GamMat);
    [draw2, lnk2, ind_red2] = fn_rmvgt_robust(M/2, mit2, kernel);
    % col 1-4 parameters, col 5:5+hp-1 eps

    draw_opt = [draw1_eps1; draw2];

    %% IS weights
    L = true;
    kernel = @(a) posterior_t_garch_whole(a, data, S, L, hyper, GamMat);
    lnk = kernel(draw_opt);

    % exp_lnd1 = 0.5*normpdf(draw_opt(:,2)).*dmvgt(draw_opt(:,1:4),mit1,false);
    % ep = duvt(draw_opt(:,5:5+hp-1), nu, hp, false);  % false: exp dnesity
    exp_lnd1 = 0.5*sum(duvt(draw_opt(:,5:5+hp-1), draw_opt(:,4), hp, false),2).*dmvgt(draw_opt(:,1:4), mit1, false, GamMat);
    exp_lnd2 = 0.5*dmvgt(draw_opt, mit2, false, GamMat);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd = log(exp_lnd);
    

    w_opt = fn_ISwgts(lnk, lnd, false);

    %% VaR and ES IS estimates
    % % h_opt = S*(1-draw_opt(:,1)) + draw_opt(:,1)*(y_T^2);
    % % y_opt = draw_opt(:,2).*sqrt(h_opt);
    h_T = volatility_t_garch(draw_opt(:,1:4), data, S);
    [y_opt, ~] = predict_t_garch(draw_opt(:,1:4), y_T, S, h_T, hp, draw_opt(:,5:5+hp-1));
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(sim,1) = IS_estim(1,1);
    ES_IS(sim,1) = IS_estim(1,2);
    
    save(['results/sv_mitisem_',num2str(hp),'_sim.mat'],'draw_opt','lnk','lnk2','VaR_IS','ES_IS');

    fprintf('MitISEM IS VAR estimate: %6.4f. \n',VaR_IS(sim,1));
    fprintf('MitISEM IS ES estimate: %6.4f. \n',ES_IS(sim,1));

end

mean_VaR_IS = mean(VaR_IS(VaR_IS<0));
mean_ES_IS = mean(ES_IS(ES_IS<0));

NSE_VaR_IS = std(VaR_IS(VaR_IS<0));              
NSE_ES_IS = std(ES_IS(ES_IS<0));                
           

fprintf('MitISEM IS VAR (mean) estimate: %6.4f. \n',mean_VaR_IS);
fprintf('MitISEM IS ES (mean) estimate: %6.4f. \n',mean_ES_IS);

fprintf('MitISEM IS NSE VaR estimate: %6.4f. \n',NSE_VaR_IS);
fprintf('MitISEM IS NSE ES estimate: %6.4f. \n',NSE_ES_IS);

cont = cont2;

model = 'sv';
s='mitisem';
gen_out

clear GamMat x
save(['results/sv_mitisem_',num2str(hp),'.mat']);