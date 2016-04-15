function [pmit, summary] = PMitISEM(kernel_init, kernel, fn_const_X, mu_init, partition, cont, GamMat)
    d = size(mu_init,2);

    S = length(partition);
 
    N = cont.mit.N;
    Hmax = cont.mit.Hmax;
    CV_tol = cont.mit.CV_tol;
    CV_old = cont.mit.CV_old;
    norm = cont.mit.norm;
     
    resampl_on = cont.resmpl_on;
    
%% Step 0: Initialization
    % get/define initial mit density shape scale and degrees of freedom
    % aka: naive proposal density
    if isa(kernel_init, 'function_handle')
        [mu, Sigma] = fn_initopt(kernel_init, mu_init);

        try
            r = fn_testSigma(Sigma);
        catch
            r = 1;
        end

        if (r == 1)
            Sigma = ones(d,d) + 2*diag(ones(d,1));
            Sigma = reshape(Sigma,1,d^2);
            display('Initial optimzation FAILED.')        
        else
            display('Initial optimzation OK.')
        end
        mit_init.mu = mu;
        mit_init.Sigma = Sigma;
        mit_init.df = 1;
        mit_init.p = 1;
        
    elseif isa(kernel_init, 'struct')
        mit_init = kernel_init;
    end
    
    % get draws and IS weights from naive  
    [theta, lnk] = fn_rmvgt_robust(N, mit_init, kernel, resampl_on);
    lnd = dmvgt(theta, mit_init, true, GamMat);
    w = fn_ISwgts(lnk, lnd, norm);
    [CV, ~] = fn_CVstop(w, CV_old, CV_tol);
    
%% Step 1: Adaptation 
    % update scale and location using IS with draws from the naive
    % fixing df 
      
    [mu_adapt, Sigma_adapt] = fn_muSigma(theta, w);
    mit_adapt.mu = mu_adapt;
    mit_adapt.Sigma = Sigma_adapt;
    mit_adapt.df = cont.mit.dfnc;
    mit_adapt.p = 1;
    
    % get draws and IS weights from adapted  
    [theta, lnk, ~] = fn_rmvgt_robust(N, mit_adapt, kernel, resampl_on);
    lnd = dmvgt(theta, mit_adapt, true, GamMat);
    w = fn_ISwgts(lnk, lnd, norm);
    [CV_new, ~] = fn_CVstop(w, CV_old, CV_tol);
    CV = [CV, CV_new];

%% Step 2 & 3: ITERATE WITH ISEM
    iter = 0;
    hstop = false;

    % Initialise pmit: split the mit into pmit
    pmit(S) = mit_adapt; % preallocate
    Sigma_adapt = reshape(mit_adapt.Sigma,d,d);
    for s = 1:S
        [s1, s2] = fn_partition_ends(partition, d, s);
        pmit(s).mu = mit_adapt.mu(s1:s2);
        pmit(s).Sigma = reshape(Sigma_adapt(s1:s2,s1:s2),1,(s2-s1+1)^2);
        pmit(s).df = mit_adapt.df;
        pmit(s).p = mit_adapt.p;
    end

    while ((iter < iter_max) && (hstop == false))       
        
        %% Step 2: DO FOR SUBSETS
        for s = 1:S
            H_s = size(pmit(s).p,2);
            hstop_s = false;
            
            % Step 2a: APPLY ISEM to the subset s to update pmit(s)
            % use the sample theta drawn from pmit (g0) and weigts w (w0)
            [s1, s2] = fn_partition_ends(partition, d, s);
%             theta_s = theta(:,s1:s2);
            d_s = s2-s1+1;
            if (s == 1)
                % ISEM for mu
                [pmit(s), summary_new] = fn_optimt(theta(:,s1:s2), pmit(s), w, cont, GamMat);
            else
                % ISEM for beta
                % X is a function of thetas from previous components      
                X = fn_const_X(theta(:,1:s1-1)); % mu = beta*X  
                r = size(X,2);
                beta = fn_beta(theta(:,s1:s2),w,X); 
                pmit(s).mu = reshape(beta,1,r*d_s); % x = theta(:,s1:s2);
                [pmit(s), summary_new] = fn_Poptimt(theta(:,s1:s2), pmit(s), w, cont, GamMat, X);

            end
            
            % Step 2b & 2c:
            while ((H_s < Hmax) && (hstop_s == false))
                H_s = H_s+1;
                
                % select the largest weights and corresponding draws               
                lnd_current = 0;
                for ss = 1:S
                    [s1, s2] = fn_partition_ends(partition, d, ss);
                    lnd_current = lnd_current + dmvgt(theta(:,s1:s2), pmit(ss), true, GamMat);
                end
                w_current = fn_ISwgts(lnk, lnd_current, norm);
                ind_w = fn_select(w_current, cont.mit.ISpc);
                theta_nc = theta(ind_w,:);
                w_nc = w_current(ind_w);

                % NEW COMPONENT
                % compute new component's mode and scale from IS weights
                mit_nc.p = cont.mit.pnc;
                mit_nc.df = cont.mit.dfnc;
                [mit_nc.mu, mit_nc.Sigma] = fn_muSigma(theta_nc(:,s), w_nc);
                
                % COBINE OLD AND NC
                % combine the old mixture mit_new and the new component mit_nc
                pmit_old = pmit(s);
                pmit(s) = fn_updateMit(pmit(s), mit_nc); 

                % UPDATE COMBINED
                % update mode, scale and df  of all mixture components
                [pmit(s), summary_new] = fn_optimt(theta(:,s), pmit(s), w, cont, GamMat);

                H_s = size(pmit(s).p,2);           
                CoV_s = ??;
                hstop_s = ??;
            end
        end
        iter = iter + 1;
        
        %% Step 3: CHECK THE UPDATED 
        % DRAW FROM UPDATED
        % get new draws from mit and evaluate new IS weights
        [theta, lnk] = fn_rmvgt_robust(N, mit_new, kernel, resampl_on);

        lnd = dmvgt(theta, mit_new, true, GamMat);
        w = fn_ISwgts(lnk, lnd, norm);

        % evaluate convergence
        CV_old = CV(size(CV,2));
        
        [CV_new, hstop_new] = fn_CVstop(w, CV_old, CV_tol);
        CV = [CV, CV_new]
        if (H > 1)
            hstop = hstop_new;
        end       
    end
    
%% Step 2: APPLY ISEM
    % optimize mixture using IS weighted EM and get draws from the new mit
    % optimize mode, scale and df
    cont.df.opt = true;
    mit_old = mit_adapt;
    [mit_new, summary_adapt] = fn_optimt(theta, mit_adapt, w, cont, GamMat);

    % get draws and log kernel evaluation
    [theta, lnk] = fn_rmvgt_robust(N, mit_new, kernel, resampl_on);
    lnd = dmvgt(theta, mit_new, true, GamMat);
    w = fn_ISwgts(lnk, lnd, norm);
	
	% stopping criteria
    [CV_new, ~] = fn_CVstop(w, CV_old, CV_tol);
    CV = [CV, CV_new];    
    
    H = length(mit_new.p);  % number of components

%% Step 3: Iterate on the number of mixture components
    % add more mixture components until convergence
    hstop = false;
    while ((H < Hmax) && (hstop == false))
        H = H+1;
        % select the largest weights and corresponding draws
        ind_w = fn_select(w, cont.mit.ISpc);
        theta_nc = theta(ind_w,:);
        w_nc = w(ind_w);

        % NEW COMPONENT
        % compute new component's mode and scale from IS weights
        mit_nc.p = cont.mit.pnc;
        mit_nc.df = cont.mit.dfnc;
        [mit_nc.mu, mit_nc.Sigma] = fn_muSigma(theta_nc, w_nc);

        % COBINE OLD AND NC
        % combine the old mixture mit_new and the new component mit_nc
        mit_old = mit_new;
        mit_new = fn_updateMit(mit_new, mit_nc); 

%%% ??? %%%        
        % DRAW FROM COMBINED
        % get draws and log kernel evaluation from the new mixture mit_new 
        [theta, lnk] = fn_rmvgt_robust(N, mit_new, kernel, resampl_on);

        lnd = dmvgt(theta, mit_new, true, GamMat);
        w = fn_ISwgts(lnk, lnd, norm);
%%%%%%%%%%%
        % UPDATE COMBINED
        % update mode, scale and df  of all mixture components
        [mit_new, summary_new] = fn_optimt(theta, mit_new, w, cont, GamMat);
        H = size(mit_new.p,2);

        % DRAW FROM UPDATED
        % get new draws from mit and evaluate new IS weights
        [theta, lnk] = fn_rmvgt_robust(N, mit_new, kernel, resampl_on);

        lnd = dmvgt(theta, mit_new, true, GamMat);
        w = fn_ISwgts(lnk, lnd, norm);

        % evaluate convergence
        CV_old = CV(size(CV,2));
        
        [CV_new, hstop_new] = fn_CVstop(w, CV_old, CV_tol);
        CV = [CV, CV_new]
        if (H > 1)
            hstop = hstop_new;
        end       
    end

%%
    summary.init = summary_adapt;
    if (Hmax > 1)
        summary.EM = summary_new;
    end
    summary.CV = CV; 
    
        
    if (CV(end) > CV(end-1))
        mit_new = mit_old;
    end
end