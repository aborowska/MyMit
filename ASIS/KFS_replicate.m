% Initialisation
% clc

clear all
close all
addpath(genpath('include/'));
addpath(genpath('ASIS/'));

v_new = ver('symbolic');
v_new = v_new.Release;
if strcmp(v_new,'(R2014b)')
    v_new = 1;
else
    v_new = 0;
end

model = 'sv';
sim = 1;
algo = 'MitISEM';
M = 10000;

% Constants
x_gam = (0:0.00001:100)' + 0.00001;
GamMat = gamma(x_gam);

% Priors <-- !!!
    % prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -2.5*log(0.025), -log(gamma(2.5))];
    prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  2.5*log(0.025), -log(gamma(2.5))];
    % prior_const = [-0.5*log(2*pi), - log(beta(20, 1.5)),  -0.5*log(2), -log(gamma(0.5))];
    logpdf_norm = @(x) prior_const(1,1) -0.5*(x.^2);
    logpdf_beta = @(x) prior_const(1,2)  + (20-1)*log(x) + (1.5-1)*log(1-x); 
    % % logpdf_gamma = @(x) prior_const(1,3) + prior_const(1,4) + (2.5-1)*log(x) - x/0.025;
    logpdf_invgamma = @(x) prior_const(1,3) + prior_const(1,4) - (2.5+1)*log(x) - 0.025./x;
    % % logpdf_chi2 = @(x) prior_const(1,3) + prior_const(1,4) -0.5*log(x) - 0.5*x;

% Parameters
mu_low = [-10, 0.95, 0.01];
mu_high = [-10, 0.95, 0.25];

d = 3;

Sigma = eye(d);
Sigma = reshape(Sigma, 1, d^2);

% Simulated data: only the centred parametrisation
T = 2000;
y_low = gen_sv(T, mu_low, 'C');
y_high = gen_sv(T, mu_high, 'C');

par_NAIS_init.b = zeros(T,1);
par_NAIS_init.C = ones(T,1); 

% Initialisation
EMitISEM_Control
cont.resmpl_on = false;
cont.mit.CV_tol = 0.1;
cont.mit.N = 10000;
 
mit_init_low.mu = mu_low;
mit_init_low.Sigma = Sigma;
mit_init_low.p = cont.mit.pnc;
mit_init_low.df = cont.mit.dfnc;
    
mit_init_high = mit_init_low;
mit_init_high.mu = mu_high;

% EMitISEM approximation
kernel_low = @(a) posterior_sv(y_low, a, par_NAIS_init, prior_const, cont.nais);
kernel_high = @(a) posterior_sv(y_high, a, par_NAIS_init, prior_const, cont.nais);

[mit_low, theta_low, x_low, w_low, lnk_low, lng_y_low, lnw_x_low, CV_low] = EMitISEM(mit_init_low, kernel_low, cont, GamMat);
cont.mit.Hmax = 2;
[mit_high, theta_high, x_high, w_high, lnk_high, lng_y_high, lnw_x_high, CV_high] = EMitISEM(mit_init_high, kernel_high, cont, GamMat);

save('replicate_mit_low.mat', 'y_low','mu_low','Sigma','cont','mit_init_low',...
    'mit_low', 'theta_low', 'x_low', 'w_low', 'lnk_low', 'lng_y_low', 'lnw_x_low', 'CV_low');

save('replicate_mit_high.mat', 'y_high','mu_high','Sigma','cont','mit_init_high',...
    'mit_high', 'theta_high', 'x_high', 'w_high', 'lnk_high', 'lng_y_high', 'lnw_x_high', 'CV_high');

% MH with the cadidate from EMitISEM
[theta_low, x_low, lnw_low, lnk_low, ~, ~, ~, accept_low] = EMit_MH(M+1000, d, kernel_low, mit_low, GamMat, true);
fprintf('(Low) MH acceptance rate: %6.4f. \n',accept_low);

save('replicate_MH_low.mat', 'y_low','mu_low','Sigma','cont','mit_init_low','M','mit_low',...
     'theta_low', 'x_low', 'lnw_low', 'lnk_low', 'accept_low');
   
[theta_high, x_high, lnw_high, lnk_high, ~, ~, ~, accept_high] = EMit_MH(M+1000, d, kernel_high, mit_high, GamMat, true);
fprintf('(High) MH acceptance rate: %6.4f. \n',accept_low);
    
save('replicate_MH_high.mat', 'y_high','mu_high','Sigma','cont','mit_init_high','M','mit_high',...
     'theta_high', 'x_high', 'lnw_high', 'lnk_high', 'accept_high');

plot_on = true;
print_on = false;

theta = theta_low;
SV_autocorr;
name = 'replicate_low_autorr.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')


theta = theta_high;
SV_autocorr;
name = 'replicate_high_autorr.png';
fig = gcf;
fig.PaperPositionMode = 'auto';
print(name,'-dpng','-r0')