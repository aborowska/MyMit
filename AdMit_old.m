function [mit, summary] = AdMit_old(kernel_init, kernel, mu_init, cont, GamMat)
    
    d = size(mu_init,2);
    m_all = [];
        
    Ns = cont.Ns;
    Np = cont.Np;
    Hmax = cont.Hmax;
    
    CV_tol = cont.CV_tol;
    
    resampl_on = cont.resampl_on;

%% Step 0: Initialization
    % compute the mode and scale of the first Student-t distribution in the
    % mixture as the mode if the log kernel function and minus the Hessian
    % of the log kernel evaluated at its mode 
    % draw a set of N points theta from this first stage conadidate
    % density with small number of degrees of freedom to allow for fat tails
    % --> fixed to one df
    if isa(kernel_init, 'function_handle')
        [mu, Sigma] = fn_initopt(kernel_init, mu_init);

        r = fn_testSigma(Sigma);
        if (r == 1)
            display('Initial optimzation FAILED.')        
        else
            display('Initial optimzation OK.')
        end

        % the first component
        mit.mu = mu;
        mit.Sigma = Sigma;
        mit.df = cont.dfnc;
        mit.p = 1;
    elseif isa(kernel_init, 'struct')
        mit = kernel_init;
    end
%% Step 1: evaluate the distribution of weights    
    % get draws and IS weights from INITIAL (aka naive)     
    [theta, lnk, ~] = fn_rmvgt_robust(Ns, mit, kernel, resampl_on);
    draws = theta; % draws from the current mixture
    lnd = dmvgt(theta, mit, true, GamMat);      
    w = fn_ISwgts(lnk, lnd, false); % false: not normalized
    CV = fn_CVstop(w);
    lnd = lnd(1:Np,:);         
    
%% Step 2a: Iterate on the number of mixture components
    % add more mixture components until convergence
    H = 1;
    hstop = false;

    while ((H < Hmax) && (hstop == false))
        H = H+1;
        display(H)
        if (cont.IS.opt) % use importance sampling
%             [mu, ~] = fn_muSigma(draws, w);
            [mu, Sigma] = fn_optIS(draws, w, cont);
%             [mu, Sigma] = fn_optIS(draws, w, cont, [], true);
        else % use standard optimization
            [~, ind] = max(w);
            mu_init = draws(ind,:); % first starting value, corresponding to the highest weight
            [mu, Sigma] = fn_muSigma(draws, w);
            
            mit_up.mu = mu;
            mit_up.Sigma = Sigma;
            mit_up.df = cont.dfnc;           
            mit_up.p = 1;
            
            % get draws and IS weights from UPDATED   
            [theta_up, lnk_up, ~] = fn_rmvgt_robust(Ns, mit_up, kernel, resampl_on);
            lnd_up = dmvgt(theta_up, mit_up, true, GamMat);
            w_up = fn_ISwgts(lnk_up, lnd_up, false); % false: not normalized
            [~, ind] = max(w_up);
            mu_up = theta_up(ind,:);    
            
            % TRY: optimise using each starting value
            fn_w = @(a) -fn_weights(a, kernel, mit, true, GamMat);
            [mu_init, Sigma_init, val_init] = fn_initopt(fn_w, mu_init);
            [mu_up, Sigma_up, val_up] = fn_initopt(fn_w, mu_up); 

            % CATCH: use IS is problems with convergence
            r = [fn_testSigma(Sigma_init), fn_testSigma(Sigma_up)];
            ind_r = (r==0);
            val = [val_init, val_up]; 
            ind_v = (isfinite(val));
            ind = (ind_r & ind_v);
            if (sum(ind) == 0)
                [mu, Sigma] = fn_optIS(draws, w, cont);
            elseif (sum(ind) == 2)
%                 [~,ind] = max(val);
                [~,ind] = min(val);
                if (ind == 1)
                    mu = mu_init;
                    Sigma = Sigma_init;
                else
                    mu = mu_up;
                    Sigma = Sigma_up;
                end
            elseif (find(ind == 1) == 1)
                mu = mu_init;
                Sigma = Sigma_init;
            else
                mu = mu_up;
                Sigma = Sigma_up;
            end
        end
        %%
        mit_h.mu = mu; 
        mit_h.df = cont.dfnc;
        mit_h.p = 1;
        
        CV_new = Inf; % new - the current best
        m_new = 0;
        
        %% loop over scaling factors in matrix Sigma 
        % if they were used: so when IS was employed
        for m = 1:size(Sigma,1)
            fprintf('m = %i\n', m);
            Sigma_h = Sigma(m,:);
            mit_h.Sigma = Sigma_h; % single new current component
            % draw from the new component
            [theta_h, lnk_h, ~] = fn_rmvgt_robust(Ns, mit_h, kernel, resampl_on);

            theta_up = [theta, theta_h];
            lnk_up = [lnk, lnk_h]; 
            
            %% form matrix used in optimasation of probabilities
            lnd_up = zeros(Np, H^2);
            % copy old evaluations (on the old components)
            for ii = 1 : (H-1)
                pos = (1:(H-1)) + (ii-1)*(H-1); % seq from 1 to H-1, form H to 2H-2, from 2H-1 to 3H-3, ...
                pos_up = (1:(H-1)) + (ii-1)*H; % seq from 1 to H-1, form 1+H to 2H-1, from 1+2H to 3H-1, ...
                % in pos_up there are H, 2H, 3H, ... columns missing -
                % will be filled with evaluations of all the draws on
                % the current component
                lnd_up(:,pos_up) = lnd(:,pos);
            end
            % evaluate all the draws (previous draws and the new draw) with the current
            % Sigma (current component)
            for ii = 1 : H
                pos = ii*(1:d);
                % on pos = H, 2H, 3H, ... - evaluations of all the draws on the current
                % on pos = H - of the draw from the first component
                % on pos = 2H - of the draw from the second component,
                % ...
                lnd_up(:,ii*H) = dmvgt(theta_up(1:Np,pos), mit_h, true, GamMat);
            end
            % evaluate the new draw on each previous component
            for ii = 1 : (H-1)
                mit_ii = struct('mu',mit.mu(ii,:),'Sigma',mit.Sigma(ii,:),'p',1,'df',cont.dfnc);
                pos = ii + (H-1)*H; % 1+(H-1)H, ..., (H-1)+(H-1)H = (H-1)(H+1) = H^2-1
                % on pos = H^2 there's already evaluation of the new
                % draw on the current copmonent
                lnd_up(:,pos) = dmvgt(theta_h(1:Np,:), mit_ii, true, GamMat);
            end
            % on pos = 1:H - evalutations of the draw from the first
            % compontent on the first, the second, ..., the new
            % component
            % ...
            % on pos = 1+(H-1)H : H^2 - evaluations of the draw from
            % the new component on the first, the second, ..., the new
            % component
            
            %% Step 2b: optimize the mixing probabilities
            % optimization of probabilities
            p_up = fn_optProb(mit.p, lnk_up(1:Np,:), lnd_up, cont);
            
            % the new mixture
            mit_up = fn_updateMit(mit, mit_h, p_up);
            
            %% draw from the new mixture - use the already available
            % computations
            comp = datasample(1:H,Ns,'Weights',p_up);
            draws_up = zeros(Ns,d);
            lnk_mix = zeros(Ns,1);
            tail = 0;
            for ii = 1:H
                nh = length(comp(comp == ii));          
                if (nh > 0)
                    pos = d*ii + (1-d:0);
                    draws_up(tail+1:tail+nh,:) = theta_up(1:nh,pos);
                    lnk_mix(tail+1:tail+nh,:) = lnk_up(1:nh,ii);
                end
                tail = tail + nh;
            end           
%%            
            % log of mixture for these draws
            lnd_mix = dmvgt(draws_up, mit_up, true, GamMat);       
            w_up = fn_ISwgts(lnk_mix, lnd_mix, false); % false: not normalized
            CV_up = fn_CVstop(w_up);
            
            % if CoV better than before
            if (CV_up < CV_new)          
                CV_new = CV_up; % update the current best
                w_new = w_up;
                p_new = p_up;
                Sigma_new = Sigma_h;
                draws_new = draws_up;
                theta_new = theta_up;
                lnd_new = lnd_up;
                lnk_new = lnk_up;
                m_new = m;
            end
            
        end
         
        CV_old = CV(length(CV));
        m_all = [m_all, m_new]
        CV = [CV, CV_new] % update the history of CV
        mit_h.Sigma = Sigma_new;
        % add the new component to the mixture
        mit = fn_updateMit(mit, mit_h, p_new);
        w = w_new;
        draws = draws_new;
        theta = theta_new;
        lnk = lnk_new;
        lnd = lnd_new;
        hstop = (abs((CV_new - CV_old)/CV_old) <= CV_tol);
        
    end
    
    summary.CV = CV;
    summary.m_all = m_all;
end 
    
