function [mit_new, theta, x, w_norm, lnk, lng_y, lnw_x, CV] = EMitISEM(y, mit_init, kernel_prior, kernel, cont, GamMat)
    T = size(y,1);
    N = cont.mit.N;
    N = 1000*ceil(N/1000);
    Hmax = cont.mit.Hmax;
    CV_tol = cont.mit.CV_tol;
    CV_old = cont.mit.CV_old;
    norm = cont.mit.norm;
 
%% Step 1: Initialisation
    [theta, ~] = fn_rmvgt_robust(N, mit_init, kernel_prior);

    if (N <= 2000)
        [theta, lnk] = fn_rmvgt_robust(N, mit_init, kernel, theta);
    else
        lnk = zeros(N,1);
        for ii = 1:(N/1000)
            ind = (1:1000) + (ii-1)*1000; 
            [theta(ind,:), lnk(ind,:)] = fn_rmvgt_robust(1000, mit_init, kernel, theta(ind,:));
        end
    end

    lnd = dmvgt(theta, mit_init, true, GamMat);
    w_norm = fn_ISwgts(lnk, lnd, norm);
    [CV, ~] = fn_CVstop(w_norm, CV_old, CV_tol);
    

%% Step 2: Adaptation
    [mu_adapt, Sigma_adapt] = fn_muSigma(theta, w_norm);
    mit_adapt.mu = mu_adapt;
    mit_adapt.Sigma = Sigma_adapt;
    mit_adapt.df = cont.mit.dfnc;
    mit_adapt.p = 1;

    [theta, ~] = fn_rmvgt_robust(N, mit_adapt, kernel_prior);
    if (N <= 2000)
        [theta, lnk] = fn_rmvgt_robust(N, mit_adapt, kernel, theta);  
    else
        lnk = zeros(N,1);
        for ii = 1:(N/1000)
            ind = (1:1000) + (ii-1)*1000; 
            [theta(ind,:), lnk(ind,:)] = fn_rmvgt_robust(1000, mit_adapt, kernel, theta(ind,:));
        end
    end
    
    lnd = dmvgt(theta, mit_adapt, true, GamMat);
    w_norm = fn_ISwgts(lnk, lnd, norm);
    [CV_new, ~] = fn_CVstop(w_norm, CV_old, CV_tol);
    CV = [CV, CV_new];

%% Step 3: ISEM
mit_old = mit_adapt;

    [mit_new, ~] = fn_optimt(theta, mit_adapt, w_norm, cont, GamMat);

    [theta, ~] = fn_rmvgt_robust(N, mit_new, kernel_prior);
    if (N <= 2000)
        [theta, lnk, ~, x, lng_y, lnw_x] = fn_rmvgt_robust(N, mit_new, kernel, theta);
    else
        x = zeros(N,T); 
        lnk = zeros(N,1);
        lng_y = zeros(N,1);
        lnw_x = zeros(N,1);
        for ii = 1:(N/1000)
            ind = (1:1000) + (ii-1)*1000; 
            [theta(ind,:), lnk(ind,:), ~, x(ind,:), lng_y(ind,:), lnw_x(ind,:)] = fn_rmvgt_robust(1000, mit_new, kernel, theta(ind,:));
        end
    end
    
    lnd = dmvgt(theta, mit_new, true, GamMat);
    w_norm = fn_ISwgts(lnk, lnd, norm);
    [CV_new, ~] = fn_CVstop(w_norm, CV_old, CV_tol);
    CV = [CV, CV_new];    

%% Step 4: Iteration on number of mixture components

    H = length(mit_new.p);  % number of components
    hstop = false;

    while ((H < Hmax) && (hstop == false))
        H = H+1;

        ind_nc = fn_select(w_norm,cont.mit.ISpc);
        theta_nc = theta(ind_nc,:);
        w_nc = w_norm(ind_nc);
        mit_nc.p = cont.mit.pnc;
        mit_nc.df = cont.mit.dfnc;
        [mit_nc.mu, mit_nc.Sigma] = fn_muSigma(theta_nc, w_nc);
        mit_old = mit_new;
        mit_new = fn_updateMit(mit_new, mit_nc); 

        [theta, ~] = fn_rmvgt_robust(N, mit_new, kernel_prior);

        if (N <= 2000)
            [theta, lnk, ~, ~, ~, ~] = fn_rmvgt_robust(N, mit_new, kernel, theta);
        else
            lnk = zeros(N,1);
            for ii = 1:(N/1000)
                ind = (1:1000) + (ii-1)*1000; 
                [theta(ind,:), lnk(ind,:), ~, ~, ~, ~] = fn_rmvgt_robust(1000, mit_new, kernel, theta(ind,:));
            end
        end
                
        lnd = dmvgt(theta, mit_new, true, GamMat);
        w_norm = fn_ISwgts(lnk, lnd, norm);

        %% UPDATE COMBINED
        [mit_new, ~] = fn_optimt(theta, mit_new, w_norm, cont, GamMat);

        % DRAW FROM UPDATED
        % get new draws from mit and evaluate new IS weights

        [theta, ~] = fn_rmvgt_robust(N, mit_new, kernel_prior);

        if (N <= 2000)
            [theta, lnk, ~, x, lng_y, lnw_x] = fn_rmvgt_robust(N, mit_new, kernel, theta);
        else
            x = zeros(N,T); 
            lnk = zeros(N,1);
            lng_y = zeros(N,1);
            lnw_x = zeros(N,1);
            for ii = 1:(N/1000)
                ind = (1:1000) + (ii-1)*1000; 
                [theta(ind,:), lnk(ind,:), ~, x(ind,:), lng_y(ind,:), lnw_x(ind,:)] = fn_rmvgt_robust(1000, mit_new, kernel, theta(ind,:));
            end
        end
        
        lnd = dmvgt(theta, mit_new, true, GamMat);
        w_norm = fn_ISwgts(lnk,lnd, false);
        
        [CV_new, hstop_new] = fn_CVstop(w_norm, CV_old, CV_tol);

        CV_old = CV(size(CV,2));
        [CV_new, hstop_new] = fn_CVstop(w_norm, CV_old, CV_tol);
        CV = [CV, CV_new];
        if (H > 1)
            hstop = hstop_new;
        end       
    end
    
    if (CV(end) > CV(end-1))
        mit_new = mit_old;
    end
end