clear all
N = 1000;
addpath('include/');
addpath('include/NAIS/');
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