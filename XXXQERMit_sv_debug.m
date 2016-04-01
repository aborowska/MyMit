clear all
N = 1000;
addpath(genpath('include/'));
addpath('results/');

prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  2.5*log(0.025), -log(gamma(2.5))];

x_gam = (0:0.00001:100)'+0.00001;
GamMat = gamma(x_gam);
 model = 'sv';
EMitISEM_Control

y = csvread('IBM_ret.csv');
y = 100*y;

T = size(y,1);
par_NAIS_init.b = zeros(T,1);
par_NAIS_init.C = ones(T,1); 

%% SVt
load('results/svt_mit1_2000.mat','mit1')

kernel_prior = @(a) prior_svt(a, prior_const); 
kernel = @(a) posterior_svt(y, a, par_NAIS_init, prior_const, cont.nais);

    [theta, ~] = fn_rmvgt_robust(N, mit1, kernel_prior); 
    theta = [mu_init(1:3),2.01];
    theta = [repmat(mu_init(1:3),N,1),(linspace(2.51,20,N))'];

[theta, lnk, ~, x, lng_y, lnw_x, x_smooth] = fn_rmvgt_robust(N, mit1, kernel, theta);
x=x';
% nu <= 2.721 ==> bad x_smooth!
nu = 2.01;
pdf_const = log(gamma((nu+1)/2)) - log(gamma(nu/2)) - 0.5*log(nu-2);
% theta_GH = theta_smooth(ii,1) + sqrt(V_smooth(ii,1))*z;
theta_GH = 0.006 + sqrt(0.05)*cont.nais.GH.z;

for ii=1:20
    Y = pdf_const - 0.5*(theta_GH + (nu+1)*log(1 + (y(ii)^2)./((nu-2).*exp(theta_GH)))); 
    w = sqrt(cont.nais.GH.h);
    w_theta_GH = w.*theta_GH;
    w_theta_GH2 = -0.5*w.*theta_GH.^2;
    X = [w, w_theta_GH, w_theta_GH2];
    XT = X';
    Y = w.*Y;
    beta = (XT*X)\(XT*Y)
end

% load('results/svt_mit2_2000.mat','mit2')
% kernel_prior = @(a) prior_svt_hl(a, prior_const); 
% kernel = @(a) posterior_svt_hl(y, a, VaR_prelim, par_NAIS_init, prior_const, cont.nais); 



%% SV
load('results/sv_mit1_2000.mat','mit1')
kernel_prior = @(a) prior_sv(a, prior_const); 
kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);


    [theta, ~] = fn_rmvgt_robust(N, mit1, kernel_prior); 
  
    theta = [repmat([0.4982    0.9992],N,1),linspace(0.00001,0.1,N)'];

[theta, lnk, ~, x, lng_y, lnw_x, x_smooth] = fn_rmvgt_robust(N, mit1, kernel, theta);
[~,ind10]=sort(lnk);
ind10 = ind10(1:10);
theta=theta(ind10,:);


load('results/sv_mit2_2000.mat','mit2','VaR_prelim')
kernel_prior = @(a) prior_sv_hl(a, prior_const); 
kernel = @(a) posterior_sv_hl(y, a, VaR_prelim, par_NAIS_init, prior_const, cont.nais); 

    [theta, ~] = fn_rmvgt_robust(N, mit2, kernel_prior); 
    
    theta = [linspace(-2,2,N),repmat(mu_init(2:3),N,1)];

[theta, lnk, ~, x, lng_y, lnw_x, x_smooth] = fn_rmvgt_robust(N, mit2, kernel, theta);




%% COMPARE DIRECT SAMPLING WITH FR_ROBUST
% 1. CORECTNESS OF RESULTS
% 2. SPEED

resampl_on = false;
addpath(genpath('include/'));
load('C:\Users\aba228\Dropbox\MPHIL\VU\MitISEM\MyMit\results\sv_mit1_10000_ibmdata.mat','cont','mit1')
mit = mit1;
kernel_prior = @(a) prior_sv(a, prior_const); 
kernel = @(a) posterior_sv(y, a, par_NAIS_init, prior_const, cont.nais);
N = 2000;
T = length(y);

%  EMit_MH:
profile on
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
[theta_MH, x_MH, lnw_MH, lnk_MH, lng_y_MH, lnw_x_MH, ~, accept] = EMit_MH(y, N, kernel_prior, kernel, mit, GamMat, false);
ind_MH  = find(~isnan(lng_y_MH));
 lnw_MH =lnw_MH(ind_MH,:);
 lnk_MH = lnk_MH(ind_MH,:); 
 lng_y_MH = lng_y_MH(ind_MH,:); 
 lnw_x_MH = lnw_x_MH(ind_MH,:);
profile off
profile viewer
p_MH = profile('info');
save results/profile_MH p_MH

% what's inside EMit_MH:
profile on
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

[theta, ~] = fn_rmvgt_robust(N, mit, kernel_prior, resampl_on);
[theta, lnk, ~, x, lng_y, lnw_x] = fn_rmvgt_robust(N, mit, kernel, resampl_on, theta);
ind = find(~isnan(lng_y));
 lnk = lnk(ind,:); 
 lng_y = lng_y(ind,:); 
 lnw_x = lnw_x(ind,:);
profile off
profile viewer
p = profile('info');
save results/profile p

%% if (N <= 2000)
% %         [theta, lnk, ~, x, lng_y, lnw_x, x_smooth] = fn_rmvgt_robust(N, mit, kernel, theta);
%     [theta, lnk, ~, x, lng_y, lnw_x] = fn_rmvgt_robust(N, mit, kernel, resampl_on, theta);
% else
%     lnk = zeros(N,1);
%     x = zeros(N,T); 
%     lng_y = zeros(N,1);
%     lnw_x = zeros(N,1);
% %         x_smooth = zeros(T,N);
%     for ii = 1:(N/1000)
%         ind = (1:1000) + (ii-1)*1000; 
% %             [theta(ind,:), lnk(ind,:), ~, x(ind,:), lng_y(ind,:), lnw_x(ind,:), x_smooth(:,ind)] = fn_rmvgt_robust(1000, mit, kernel, theta(ind,:));      
%         [theta(ind,:), lnk(ind,:), ~, x(ind,:), lng_y(ind,:), lnw_x(ind,:)] = fn_rmvgt_robust(1000, mit, kernel, resampl_on, theta(ind,:));      
%     end
% end
    
    
%% direct sampling
profile on
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);     

theta_direct = rmvgt2(N, mit1.mu, mit1.Sigma, mit1.df, mit1.p); 
lnk_direct = zeros(N,1);
x_direct = zeros(N,T);
lng_y_direct = zeros(N,1);
lnw_x_direct = zeros(N,1);
for ii = 1:N
    [par_NAIS, x_smooth]= NAIS_param(par_NAIS_init, y, theta_direct(ii,1:3), cont.nais);
    [x_direct(ii,:), lng_y_direct(ii,1), lnw_x_direct(ii,1)] = NAIS_loglik(y, theta_direct(ii,1:3), par_NAIS, cont.nais);
    prior = prior_sv(theta_direct(ii,1:3),prior_const);
    lnk_direct(ii,1) = lng_y_direct(ii,1) + lnw_x_direct(ii,1) + prior;
end
ind_direct  = find(~isnan(lng_y_direct));
 lnk_direct = lnk_direct(ind_direct,:); 
 lng_y_direct = lng_y_direct(ind_direct,:); 
 lnw_x_direct = lnw_x_direct(ind_direct,:);
profile off
profile viewer
p_direct = profile('info');
save results/profile_direct  p_direct

profview(0,p_MH)
profview(0,p )
profview(0,p_direct)