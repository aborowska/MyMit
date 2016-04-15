clear all
addpath('include/');
addpath('include/NAIS/');
addpath('results/');

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);
% prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -2.5*log(0.025), -log(gamma(2.5))];
prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  2.5*log(0.025), -log(gamma(2.5))];
% prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -0.5*log(2), -log(gamma(0.5))];

logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
% logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
logpdf_invgamma = @(x) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(x) - 0.025./x;
% logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) -0.5*log(x) - 0.5*x;

model = 'sv';
y = csvread('IBM_ret.csv');
y = 100*y;

T = size(y,1);
 
EMitISEM_Control

cont.mit.dfnc = 5;
N = cont.mit.N;
Hmax = cont.mit.Hmax;
CV_tol = cont.mit.CV_tol;
CV_old = cont.mit.CV_old;
norm = cont.mit.norm;
    
par_NAIS_init.b = zeros(T,1);
par_NAIS_init.C = ones(T,1); 

%% Initialisation
load('SML_ibm.mat', 'par_SV_opt', 'V_SV_corr_opt') 
d = size(par_SV_opt,2);
 
Sigma = V_SV_corr_opt;
Sigma = reshape(Sigma,1,d^2);
mit_init.mu = par_SV_opt;
mit_init.Sigma = Sigma;
mit_init.p = 1;
mit_init.df = 5;

kernel = @(a) prior_sv(a, true, prior_const); 
[theta, prior] = fn_rmvgt_robust(N, mit_init, kernel);

kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);
[theta, lnk, ~, x, lng_y, lnw_x, par_NAIS] = fn_rmvgt_robust(N, mit_init, kernel, theta);

%% fn_par_zero;
% 
% for ii = 1:N
%     if (mod(ii,500)==0)
%         fprintf('iter = %i\n',ii)
%     end
%     par_NAIS_iter = NAIS_param(par_NAIS_init, y, theta(ii,:), cont.nais);
%     par_NAIS.b(:,ii) = par_NAIS_iter.b;
%     par_NAIS.C(:,ii) = par_NAIS_iter.C;
% end
% 
% b = par_NAIS.b;
% C = par_NAIS.C;
% a = 0.5*(log(abs(C)) - log(2*pi) - (b.^2)./C);
% y_star = b./C;
% 
% for ii = 1:N
%     if (mod(ii,500)==0)
%         fprintf('iter = %i\n',ii)
%     end
%     par_NAIS_iter.b = b(:,ii);
%     par_NAIS_iter.C = C(:,ii);
%     par_KFS = IS_Model(par_NAIS_iter, theta(ii,:));
%     [~, ~, v(:,ii) , F_inv(:,ii) , eps_smooth(:,ii) , K(:,ii) , L(:,ii)] = KFS(y_star(:,ii),par_KFS); 
%     lng_y(ii,1) = -0.5*(T*log(2*pi) + sum((v(:,ii).^2).*F_inv(:,ii),1) - sum(log(abs(F_inv(:,ii))),1));
%     x_sim(:,ii) = IS_sim(y_star(:,ii), eps_smooth(:,ii) , v(:,ii) , F_inv(:,ii) , K(:,ii) , L(:,ii) , par_KFS);
%     lnp(ii,1) = sum( -0.5*(log(2*pi) +  x_sim(:,ii)  + (y.^2)./exp(x_sim(:,ii))) );
%     lng(ii,1) = sum( a(:,ii) + b(:,ii).*x_sim(:,ii) - 0.5*C(:,ii).*x_sim(:,ii).^2 );
% end
%%
lnd = dmvgt(theta, mit_init, true, GamMat);
% lnw = lng_y + lnp - lng + prior - lnd;
lnw = lng_y + lnw_x + prior - lnd;
w = exp(lnw - max(lnw));

w_norm = w/sum(w);
[CV, ~] = fn_CVstop(w_norm, CV_old, CV_tol);


%% Adaptation
[mu_adapt, Sigma_adapt] = fn_muSigma(theta, w);
mit_adapt.mu = mu_adapt;
mit_adapt.Sigma = Sigma_adapt;
mit_adapt.df = cont.mit.dfnc;
mit_adapt.p = 1;

kernel = @(a) prior_sv(a, true, prior_const); 
[theta, prior] = fn_rmvgt_robust(N, mit_adapt, kernel);

kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);
[theta, lnk, ~, x, lng_y, lnw_x, par_NAIS] = fn_rmvgt_robust(N, mit_adapt, kernel, theta);

%% fn_par_zero;
% 
% for ii = 1:N
%     if (mod(ii,500)==0)
%         fprintf('iter = %i\n',ii)
%     end
%     par_NAIS_iter = NAIS_param(par_NAIS_init, y, theta(ii,:), cont.nais);
%     par_NAIS.b(:,ii) = par_NAIS_iter.b;
%     par_NAIS.C(:,ii) = par_NAIS_iter.C;
% end
% 
% b = par_NAIS.b;
% C = par_NAIS.C;
% a = 0.5*(log(abs(C)) - log(2*pi) - (b.^2)./C);
% y_star = b./C;
% 
% for ii = 1:N
%     if (mod(ii,500)==0)
%         fprintf('iter = %i\n',ii)
%     end
%     par_NAIS_iter.b = b(:,ii);
%     par_NAIS_iter.C = C(:,ii);
%     par_KFS = IS_Model(par_NAIS_iter, theta(ii,:));
%     [~, ~, v, F_inv, eps_smooth, K, L] = KFS(y_star(:,ii),par_KFS); 
%     lng_y(ii,1) = -0.5*(T*log(2*pi) + sum((v.^2).*F_inv,1) - sum(log(abs(F_inv)),1));
%     x_sim(:,ii) = IS_sim(y_star(:,ii), eps_smooth, v, F_inv, K, L, par_KFS);
%     lnp(ii,1) = sum( -0.5*(log(2*pi) +  x_sim(:,ii)  + (y.^2)./exp(x_sim(:,ii))) );
%     lng(ii,1) = sum( a(:,ii) + b(:,ii).*x_sim(:,ii) - 0.5*C(:,ii).*x_sim(:,ii).^2 );
% end
%%
lnd = dmvgt(theta, mit_adapt, true, GamMat);
% lnw = lng_y + lnp - lng + prior - lnd;
lnw = lng_y + lnw_x + prior - lnd;
w = exp(lnw - max(lnw));

w_norm = w/sum(w);
[CV_new, ~] = fn_CVstop(w_norm, CV_old, CV_tol);
CV = [CV, CV_new];


%% ISEM
[mit_new, summary_new] = fn_optimt(theta, mit_adapt, w_norm, cont, GamMat);

kernel = @(a) prior_sv(a, true, prior_const); 
[theta, prior] = fn_rmvgt_robust(N, mit_new, kernel);
 
kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);
[theta, lnk, ~, x, lng_y, lnw_x, par_NAIS] = fn_rmvgt_robust(N, mit_new, kernel, theta);

%% fn_par_zero;
%
% for ii = 1:N
%     if (mod(ii,500)==0)
%         fprintf('iter = %i\n',ii)
%     end
%     par_NAIS_iter = NAIS_param(par_NAIS_init, y, theta(ii,:), cont.nais);
%     par_NAIS.b(:,ii) = par_NAIS_iter.b;
%     par_NAIS.C(:,ii) = par_NAIS_iter.C;
% end
% 
% b = par_NAIS.b;
% C = par_NAIS.C;
% a = 0.5*(log(abs(C)) - log(2*pi) - (b.^2)./C);
% y_star = b./C;
% 
% for ii = 1:N
%     if (mod(ii,500)==0)
%         fprintf('iter = %i\n',ii)
%     end
%     par_NAIS_iter.b = b(:,ii);
%     par_NAIS_iter.C = C(:,ii);
%     par_KFS = IS_Model(par_NAIS_iter, theta(ii,:));
%     [~, ~, v, F_inv, eps_smooth, K, L] = KFS(y_star(:,ii),par_KFS); 
%     lng_y(ii,1) = -0.5*(T*log(2*pi) + sum((v.^2).*F_inv,1) - sum(log(abs(F_inv)),1));
%     x_sim(:,ii) = IS_sim(y_star(:,ii), eps_smooth, v, F_inv, K, L, par_KFS);
%     lnp(ii,1) = sum( -0.5*(log(2*pi) +  x_sim(:,ii)  + (y.^2)./exp(x_sim(:,ii))) );
%     lng(ii,1) = sum( a(:,ii) + b(:,ii).*x_sim(:,ii) - 0.5*C(:,ii).*x_sim(:,ii).^2 );
% end
%%
lnd = dmvgt(theta, mit_new, true, GamMat);
% lnw = lng_y + lnp - lng + prior - lnd;
lnw = lng_y + lnw_x + prior - lnd;

w = exp(lnw - max(lnw));

w_norm = w/sum(w);

% [CV, ~] = fn_CVstop(w, CV_old, CV_tol);
[CV_new, ~] = fn_CVstop(w_norm, CV_old, CV_tol);
CV = [CV, CV_new];
H = length(mit_new.p);  % number of components
 
%% New components
hstop = false;
while ((H < Hmax) && (hstop == false))
    H = H+1;

    ind_nc = fn_select(w_norm,cont.mit.ISpc);
    theta_nc = theta(ind_nc,:);
    w_nc = w_norm(ind_nc);
    mit_nc.p = cont.mit.pnc;
    mit_nc.df = cont.mit.dfnc;
    [mit_nc.mu, mit_nc.Sigma] = fn_muSigma(theta_nc, w_nc);
    mit_new = fn_updateMit(mit_new, mit_nc); 

    kernel = @(a) prior_sv(a, true, prior_const); 
    [theta, prior] = fn_rmvgt_robust(N, mit_new, kernel);

    kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);
    [theta, lnk, ~, x, lng_y, lnw_x, par_NAIS] = fn_rmvgt_robust(N, mit_new, kernel, theta);

%%     fn_par_zero;
%     
%     for ii = 1:N
%         if (mod(ii,500)==0)
%             fprintf('iter = %i\n',ii)
%         end
%         par_NAIS_iter = NAIS_param(par_NAIS_init, y, theta(ii,:), cont.nais);
%         par_NAIS.b(:,ii) = par_NAIS_iter.b;
%         par_NAIS.C(:,ii) = par_NAIS_iter.C;
%     end
% 
%     b = par_NAIS.b;
%     C = par_NAIS.C;
%     a = 0.5*(log(abs(C)) - log(2*pi) - (b.^2)./C);
%     y_star = b./C;
% 
%     for ii = 1:N
%         if (mod(ii,500)==0)
%             fprintf('iter = %i\n',ii)
%         end
%         par_NAIS_iter.b = b(:,ii);
%         par_NAIS_iter.C = C(:,ii);
%         par_KFS = IS_Model(par_NAIS_iter, theta(ii,:));
%         [~, ~, v, F_inv, eps_smooth, K, L] = KFS(y_star(:,ii),par_KFS); 
%         lng_y(ii,1) = -0.5*(T*log(2*pi) + sum((v.^2).*F_inv,1) - sum(log(abs(F_inv)),1));
%         x_sim(:,ii) = IS_sim(y_star(:,ii), eps_smooth, v, F_inv, K, L, par_KFS);
%         lnp(ii,1) = sum( -0.5*(log(2*pi) +  x_sim(:,ii)  + (y.^2)./exp(x_sim(:,ii))) );
%         lng(ii,1) = sum( a(:,ii) + b(:,ii).*x_sim(:,ii) - 0.5*C(:,ii).*x_sim(:,ii).^2 );
%     end
%%
    lnd = dmvgt(theta, mit_new, true, GamMat);
%     lnw = lng_y + lnp - lng + prior - lnd;
    lnw = lng_y + lnw_x + prior - lnd;
    w = exp(lnw - max(lnw));

    w_norm = w/sum(w);

    %% UPDATE COMBINED
    [mit_new, summary_new] = fn_optimt(theta, mit_new, w_norm, cont, GamMat);

    % DRAW FROM UPDATED
    % get new draws from mit and evaluate new IS weights

    kernel = @(a) prior_sv(a, true, prior_const); 
    [theta, prior] = fn_rmvgt_robust(N, mit_new, kernel);

     kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);
    [theta, lnk, ~, x, lng_y, lnw_x, par_NAIS] = fn_rmvgt_robust(N, mit_new, kernel, theta);
   
%%     fn_par_zero;
% 
%     for ii = 1:N
%         if (mod(ii,500)==0)
%             fprintf('iter = %i\n',ii)
%         end
%         par_NAIS_iter = NAIS_param(par_NAIS_init, y, theta(ii,:), cont.nais);
%         par_NAIS.b(:,ii) = par_NAIS_iter.b;
%         par_NAIS.C(:,ii) = par_NAIS_iter.C;
%     end
% 
%     b = par_NAIS.b;
%     C = par_NAIS.C;
%     a = 0.5*(log(abs(C)) - log(2*pi) - (b.^2)./C);
%     y_star = b./C;
% 
%     for ii = 1:N
%         if (mod(ii,500)==0)
%             fprintf('iter = %i\n',ii)
%         end
%         par_NAIS_iter.b = b(:,ii);
%         par_NAIS_iter.C = C(:,ii);
%         par_KFS = IS_Model(par_NAIS_iter, theta(ii,:));
%         [~, ~, v, F_inv, eps_smooth, K, L] = KFS(y_star(:,ii),par_KFS); 
%         lng_y(ii,1) = -0.5*(T*log(2*pi) + sum((v.^2).*F_inv,1) - sum(log(abs(F_inv)),1));
%         x_sim(:,ii) = IS_sim(y_star(:,ii), eps_smooth, v, F_inv, K, L, par_KFS);
%         lnp(ii,1) = sum( -0.5*(log(2*pi) +  x_sim(:,ii)  + (y.^2)./exp(x_sim(:,ii))) );
%         lng(ii,1) = sum( a(:,ii) + b(:,ii).*x_sim(:,ii) - 0.5*C(:,ii).*x_sim(:,ii).^2 );
%     end
%%
    lnd = dmvgt(theta, mit_new, true, GamMat);
%     lnw = lng_y + lnp - lng + prior - lnd;
    lnw = lng_y + lnw_x + prior - lnd;
    w = exp(lnw - max(lnw));

    w_norm = w/sum(w);

    [CV_new, hstop_new] = fn_CVstop(w_norm, CV_old, CV_tol);

    CV_old = CV(size(CV,2));
    [CV_new, hstop_new] = fn_CVstop(w, CV_old, CV_tol);
    CV = [CV, CV_new];
    if (H > 1)
        hstop = hstop_new;
    end       
end 
 