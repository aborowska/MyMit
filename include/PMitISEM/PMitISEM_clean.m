clear all
close all

%% Initialization
s = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s); 

addpath(genpath('include/'));
plot_on = false;
save_on = false;

x_gam = (0:0.00001:50)' + 0.00001; 
GamMat = gamma(x_gam);

% Artificial, white noise data 
T = 10000;
y = randn(T,1); 
y = y - mean(y);

% sigma is the VARIANCE of the error term, i.e. y_t ~ NID(0, sigma)
sigma_init = 0.9;

% Control parameters for  MitISEM 
MitISEM_Control
cont.mit.dfnc = 5;
cont.mit.N = 10000;
cont.resmpl_on = false;
% hyper parameters for the prior for sigma2(inv. gamma)
a = 1; % if a == 0, then the flat prior is used; if a == 1, then the conjugate prior (inv. gamma)
b = 1; 

% Construct the approximation to the parameter posterior
kernel_init = @(x) - posterior_debug(x, y, a, b, true);
kernel = @(x) posterior_debug(x, y, a, b, true);
[mit1, summary1] = MitISEM_new(kernel_init, kernel, sigma_init, cont, GamMat);

B = 1;
VaR_prelim = zeros(B,1);
VaR_IS = zeros(B,1);

% Metropolis-Hastings for the parameters
M = 100000;
BurnIn = 1000;

H = 10; % forecast horizon
d = H+1; % dimension of theta
% partition = 1:H+1;
partition = [1,2:H+1];
S = length(partition); % number of subsets    


model = ['WN',sprintf('%d',partition(1,2)-partition(1,1))];
algo = 'MitISEM';
p_bar = 0.01;


for bb = 1:B    
    kernel = @(x) posterior_debug(x, y, a, b, true);
    [sigma1, accept ] = Mit_MH(M+BurnIn, kernel, mit1, GamMat);
    fprintf('(%s) MH acceptance rate: %4.2f. \n', model, accept);
    sigma1 = sigma1(1001:M+BurnIn);

    %% Future disturbances
    eps_H = randn(M,H); % --> if future disturbances drawn from the target then their weights are 1

    % the returns coresponding to disturbances: y = sqrt(sigma)*eps
    y_H = bsxfun(@times,eps_H,sqrt(sigma1)); 
    % preliminary VaR
    [PL, ind] = sort(fn_PL(y_H));
    VaR_prelim(bb,1) = PL(p_bar*M);  
    ES_prelim = mean(PL(1:p_bar*M));    
    fprintf('Preliminary 100*%4.2f%% VaR estimate: %6.4f (%s, %s). \n', p_bar, VaR_prelim(bb,1), model, algo);
    % theoretical VaR (with the FLAT prior) (or, the likelihood)
    VaR_true = norminv(p_bar,0,sqrt(H));
    fprintf('Theoretical 100*%4.2f%% VaR (under flat prior): %6.4f.\n', p_bar,VaR_true);

    % High loss draws = the target of the truncated H-days-ahead return distibution
    draw = [sigma1, eps_H];
    draw_hl = draw(ind,:);
    draw_hl = draw_hl(PL<=VaR_prelim(bb,1),:);  

    % WEIGHTS to initialise PMitISEM
    % future disturbances are generated from the target thus have weights 1
    % log kernel evaluation - only for the parameter draws
    kernel = @(x) posterior_debug(x, y, a, b, true);
    lnk = kernel(draw(:,1));
    lnk_hl = kernel(draw_hl(:,1)); 

    % log candidate evaluation
    lnd = dmvgt(draw(:,1), mit1, true, GamMat);
    lnd_hl = dmvgt(draw_hl(:,1), mit1, true, GamMat);

    % log weights
    w = lnk - lnd;
    % w = exp(w - max(w));

    w_hl = lnk_hl - lnd_hl;
    w_hl = exp(w_hl - max(w_hl));

    %% Use draw_hl and w_hl to start the Partial MitISEM

    % patial mixture - based on multidimensional stuctures with 4 mit fields
    % mit_struct = struct('mu',[],'Sigma',[],'df',[],'p',[]);
    % pmit_struct = struct('mu',cell(1,S),'Sigma',cell(1,S),'df',cell(1,S),'p',cell(1,S));

    % TOY PMIT - without interdependance between sets
    pmit = struct('mu',cell(1,S),'Sigma',cell(1,S),'df',cell(1,S),'p',cell(1,S));
    for s = 1:S
        [s1, s2] = fn_partition_ends(partition, d, s);
        theta_s = draw_hl(:,s1:s2);
        [mu_adapt, Sigma_adapt] = fn_muSigma(theta_s, w_hl);

    %     [mu_adapt, Sigma_adapt] = fnr_muSigma(draw_hl(:,s), w_hl);
        mit_adapt.mu = mu_adapt;
        mit_adapt.Sigma = Sigma_adapt;
        mit_adapt.df = cont.mit.dfnc;
        mit_adapt.p = 1;

    %     pmit(s) = fn_optimt(draw_hl(:,s), mit_adapt, w_hl, cont, GamMat);   
        pmit(s) = fn_optimt(theta_s, mit_adapt, w_hl, cont, GamMat);
    end

    % draw from pmit
    draw_pmit = zeros(M,d);

    % kernel = @(x) posterior_debug(x, y, a, b, true);
    kernel = @(x) posterior_debug_hl(x, y, a, b, Inf, true); 
    [s1, s2] = fn_partition_ends(partition, d, 1);
%     [draw_pmit(:,s1:s2), lnk] = fn_rmvgt_robust(M, pmit(1), kernel, false);
    [draw_pmit(:,s1:s2), ~] = fn_rmvgt_robust(M, pmit(1), kernel, false);
%     lnd = dmvgt(draw_pmit(:,s1:s2), pmit(1), true, GamMat);
    kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
    for ss = 2:S
        [s1, s2] = fn_partition_ends(partition, d, ss);
%         [draw_pmit(:,s1:s2), lnk_s] = fn_rmvgt_robust(M, pmit(ss), kernel, false);
        [draw_pmit(:,s1:s2), ~] = fn_rmvgt_robust(M, pmit(ss), kernel, false);
%         lnk = lnk + lnk_s;
%         lnd = lnd + dmvgt(draw_pmit(:,s1:s2), pmit(ss), true, GamMat);
    end   
%     w_pmit = lnk - lnd;

    % y_H = bsxfun(@times,eps_H,sqrt(sigma1)); 
    y_pmit = bsxfun(@times,draw_pmit(:,2:d),sqrt(draw_pmit(:,1)));  % future returns 
    % y_pmit = sort(fn_PL(return_pmit));
    % the returns coresponding to disturbances
    [PL_pmit, ~] = sort(fn_PL(y_pmit));

    y_opt = [y_H; y_pmit];
    PL_opt =  sort(fn_PL(y_opt));
    draw_opt = [draw; draw_pmit];

    kernel = @(x) posterior_debug(x, y, a, b, true);
    lnk_opt = kernel(draw_opt(:,1));

%     kernel = @(x) posterior_debug_hl(x, y, a, b, Inf, true); 
%    [draw_pmit(:,s1:s2), lnk_s] = fn_rmvgt_robust(M, pmit(ss), kernel, false);
%     lnk_opt = kernel(draw_opt(:,s1:s2));

    kernel = @(aa) - 0.5*(log(2*pi) + aa.^2);
    eps_pdf = sum(kernel(draw_opt(:,2:d)),2);
    lnk_opt = lnk_opt + eps_pdf;

    % optimal weights
    [s1, s2] = fn_partition_ends(partition, d, 1);
    exp_lnd1 = 0.5*exp(eps_pdf + dmvgt(draw_opt(:,s1:s2), mit1, true, GamMat));
    exp_lnd2 = dmvgt(draw_opt(:,s1:s2),pmit(1), true, GamMat);
    for ss = 2:S
        [s1,s2] = fn_partition_ends(partition, d, ss);
        exp_lnd2 = exp_lnd2 + dmvgt(draw_opt(:,s1:s2),pmit(ss), true, GamMat);
    end
    exp_lnd2 = 0.5*exp(exp_lnd2);
    exp_lnd = exp_lnd1 + exp_lnd2;
    lnd_opt = log(exp_lnd);
    w_opt = fn_ISwgts(lnk_opt, lnd_opt, false);

    % IS VaR estimation
    dens = struct('y',y_opt,'w',w_opt,'p_bar',p_bar);
    IS_estim = fn_PL(dens, 1);
    VaR_IS(bb,1) = IS_estim(1,1);
    ES_IS = IS_estim(1,2);   
end

if plot_on
    figure(6) 
    set(gcf,'units','normalized','outerposition',[0 0 0.5 0.5]);   
    set(gcf,'defaulttextinterpreter','latex');
    boxplot([VaR_prelim,VaR_IS],'labels',{'VaR prelim','VaR Toy PMit'})
    plotTickLatex2D;
    name = ['figures/PMitISEM/',model,'_', num2str(p_bar),'_VaR_box_B',num2str(B),'.png'];
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(name,'-dpng','-r0')
end

if plot_on
    figure(2)
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.6 0.6]);
    set(gcf,'defaulttextinterpreter','latex');
    
    subplot(1,2,1)
    hold on
    plot(PL)
    plot(sort(PL_pmit),'r')
    plot(ones(M,1)*VaR_prelim(bb,1),'c')
    plot(ones(M,1)*VaR_IS(bb,1),'m')
    hold off
    legend('Direct returns','Toy PMit returns','VaR prelim', 'VaR IS')
    plotTickLatex2D
    
    subplot(1,2,2)
    hold on 
    plot(PL_opt)
    pos =  max(find(PL_opt<=VaR_prelim(bb,1)));
    scatter(pos, VaR_prelim(bb,1),'MarkerEdgeColor','green','MarkerFaceColor','green')   
    pos =  max(find(PL_opt<=VaR_IS(bb,1)));
    scatter(pos, VaR_IS(bb,1),'MarkerEdgeColor','red','MarkerFaceColor','red')      
    hold off
    legend('Combined returns','VaR prelim','VaR IS')
    plotTickLatex2D;
    
    name = ['figures/PMitISEM/',model,'_', num2str(p_bar),'_PL',num2str(B),'.png'];
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(name,'-dpng','-r0')
end

if save_on
    name = ['figures/PMitISEM/',model,'_', num2str(p_bar),'_VaR_results',num2str(B),'.mat'];
    save(name,'VaR_prelim','VaR_IS')
end
% pmit_trial = pmit;





%% REAL PMIT
partition = [1,3:H+1];
S = length(partition);
d = H+1;

fn_const_X = @(a) WN_const_X(a);


kernel = @(x) posterior_debug_hl(x, y, a, b, VaR_prelim(bb,1), true); 
CV_old = cont.mit.CV_old;
CV_tol = cont.mit.CV_tol;

% Initialization - just rejection sampling
draw0 = draw_hl;
w0 = w_hl;
lnk0 = kernel(draw0);

% Adapted patial mixture - to have an input to ISEM
pmit_adapt = struct('mu',cell(1,S),'Sigma',cell(1,S),'df',cell(1,S),'p',cell(1,S));
for s = 1:S    
    [s1, s2] = fn_partition_ends(partition, d, s);
    theta_s = draw0(:,s1:s2);
    if (s==1)
        [pmit_adapt(s).mu, pmit_adapt(s).Sigma] = fn_muSigma(theta_s, w0);
    else
        X = fn_const_X(draw0(:,1:s1-1));
        %         r = size(X,2);
        %        beta = fn_beta(theta_s,w0,X); 
        %        pmit(s).mu = reshape(beta,1,r*d_s); % x = theta(:,s1:s2);
        [beta, Sigma] = fn_beta(theta_s,w0,X); 
        pmit_adapt(s).mu = beta; % x = theta(:,s1:s2);
        pmit_adapt(s).Sigma = Sigma;
    end   
    pmit_adapt(s).df = cont.mit.dfnc;
    pmit_adapt(s).p = 1;
end
[CV, ~] = fn_CVstop(w0, CV_old, CV_tol);




%% ISEM optimized partial mixture (for the time being - one component, no iteration)
pmit = pmit_adapt;
for s = 1:S    
    [s1, s2] = fn_partition_ends(partition, d, s);
    d_s = s2-s1+1;
    theta_s = draw0(:,s1:s2);
    if (s==1)
        pmit(s) = fn_optimt(theta_s, pmit(s), w0, cont, GamMat);
    else
       X = fn_const_X(draw0(:,1:s1-1));
%        r = size(X,2);        
%        beta = fn_beta(theta_s,w0,X); 
%        pmit(s).mu = reshape(beta,1,r*d_s); % x = theta(:,s1:s2);
       [beta, Sigma] = fn_beta(theta_s,w0,X); 
       pmit(s).mu = beta; % x = theta(:,s1:s2);
       pmit(s).Sigma = Sigma;
       pmit(s) = fn_Poptimt(theta_s, pmit(s), w0, cont, GamMat, X);
    end
%     lnd_curr = fn_dpmit(draw0, pmit, partition, true, GamMat);
    % theta = draw0; L = true;
    lnd_curr = fn_dpmit(draw0, pmit, partition, fn_const_X, true, GamMat);
    
    w_curr = fn_ISwgts(lnk0, lnd_curr, false); % we keep lnk0 fixed, only candidate evaluation changes
    [CV_new, ~] = fn_CVstop(w_curr, CV_old, CV_tol);
    CV = [CV, CV_new]; 
    % INCREASE NECESSARY: (??) 
    % with w0 I only evaluated sigmas (parameter posterior vs mit1)
    % with w_curr I evaluate d parameters: sigma + H epsilons
end

pmit_one = pmit;

%% Adding more components (for the time being - without step 3)
% ISEM optimized partial mixture (with iterations)

CV = cell(S,1);
% CV(:) = {CV_old};
[CV0, ~] = fn_CVstop(w0, CV_old, CV_tol);
CV(:) = {CV0};
Hmax = cont.mit.Hmax;
pmit = pmit_adapt;
for s = 1:S 
    hstop = false;
    H_s = 1;  
    
    [s1, s2] = fn_partition_ends(partition, d, s);
    d_s = s2-s1+1;
    theta_s = draw0(:,s1:s2);
    if (s==1)
        pmit(s) = fn_optimt(theta_s, pmit(s), w0, cont, GamMat);
    else
       X = fn_const_X(draw0(:,1:s1-1));
       [beta, Sigma] = fn_beta(theta_s,w0,X); 
       pmit(s).mu = beta; % x = theta(:,s1:s2);
       pmit(s).Sigma = Sigma;
       pmit(s) = fn_Poptimt(theta_s, pmit(s), w0, cont, GamMat, X);
    end

    lnd_curr = fn_dpmit(draw0, pmit, partition, fn_const_X, true, GamMat);
    
    w_curr = fn_ISwgts(lnk0, lnd_curr, false); % we keep lnk0 fixed, only candidate evaluation changes
    CV_old =  CV{s};
    [CV_new, ~] = fn_CVstop(w_curr, CV_old, CV_tol);
    CV(s) = {[CV{s}, CV_new]}; 

    % Step 2b & 2c:
    while ((H_s < Hmax) && (hstop == false))
        H_s = H_s+1;
        ind_w = fn_select(w_curr, cont.mit.ISpc);
        theta_nc = draw0(ind_w,:);
        w_nc = w_curr(ind_w);
        
        theta_nc_s = theta_nc(:,s1:s2);

        % NEW COMPONENT
        % compute new component's mode and scale from IS weights
        mit_nc.p = cont.mit.pnc;
        mit_nc.df = cont.mit.dfnc;
        if (s==1)
            [mit_nc.mu, mit_nc.Sigma] = fn_muSigma(theta_nc_s, w_nc);
        else
            X = fn_const_X(theta_nc(:,1:s1-1));
            [beta, Sigma] = fn_beta(theta_nc_s,w_nc,X); 
            mit_nc.mu = beta; 
            mit_nc.Sigma = Sigma;
        end              
                        
        % COBINE OLD AND NC
        % combine the old mixture mit_new and the new component mit_nc
        pmit_old = pmit(s);
        pmit(s) = fn_updateMit(pmit(s), mit_nc); 
        
        % UPDATE COMBINED
        % update mode, scale and df  of all mixture components
        if (s==1)
            pmit(s) = fn_optimt(theta_s, pmit(s), w0, cont, GamMat); % w_curr only to locate the new component; w0 to run ISEM
        else
           X = fn_const_X(draw0(:,1:s1-1));
%            [beta, Sigma] = fn_beta(theta_s,w0,X); 
%            pmit(s).mu = beta; 
%            pmit(s).Sigma = Sigma;
           pmit(s) = fn_Poptimt(theta_s, pmit(s), w0, cont, GamMat, X);
        end

        lnd_curr = fn_dpmit(draw0, pmit, partition, fn_const_X, true, GamMat);  % for all draws!

        w_curr = fn_ISwgts(lnk0, lnd_curr, false); % we keep lnk0 fixed, only candidate evaluation changes
        CV_old = CV{s}; 
        CV_old = CV_old(:,end);
        [CV_new, hstop_new] = fn_CVstop(w_curr, CV_old, CV_tol);
        CV(s) = {[CV{s}, CV_new]}; 
%         if (H > 1)
            hstop = hstop_new;
%         end   
    end
    
%     CV_old = CV{s}; 
%     if (CV_old(:,end) > CV_old(:,end-1))
%         pmit(s) = pmit_old;
%     end
    
end


if save_on
    name = ['figures/PMitISEM/',model,'_', num2str(p_bar),'_pmits.mat'];
    save(name,'pmit_adapt','pmit_one','pmit','CV')
end


%% ADD STEP 3: sample from pmit
kernel = @(x) posterior_debug_hl(x, y, a, b, VaR_prelim(bb,1), true); 
[draw_pmit ,lnk_pmit] = fn_p_rmvgt(size(draw0,1), pmit, d, partition, kernel, fn_const_X);  
lnd_pmit = fn_dpmit(draw_pmit, pmit, partition, fn_const_X, true, GamMat);
w_pmit = fn_ISwgts(lnk_pmit, lnd_pmit, false); 
[CV_pmit, ~] = fn_CVstop(w_pmit, CV_old, CV_tol);

draw0 = draw_pmit;
w0 = w_pmit;
lnk0 = lnk_pmit;
